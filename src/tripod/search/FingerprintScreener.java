package tripod.search;

/***********************************************************************
			 PUBLIC DOMAIN NOTICE
		     NIH Chemical Genomics Center
         National Center for Advancing Translational Sciences

This software/database is a "United States Government Work" under the
terms of the United States Copyright Act.  It was written as part of
the author's official duties as United States Government employee and
thus cannot be copyrighted.  This software/database is freely
available to the public for use. The NIH Chemical Genomics Center
(NCGC) and the U.S. Government have not placed any restriction on its
use or reproduction. 

Although all reasonable efforts have been taken to ensure the accuracy
and reliability of the software and data, the NCGC and the U.S.
Government do not and cannot warrant the performance or results that
may be obtained by using this software or data. The NCGC and the U.S.
Government disclaim all warranties, express or implied, including
warranties of performance, merchantability or fitness for any
particular purpose.

Please cite the authors in any work or product based on this material.

************************************************************************/

import java.io.PrintStream;

import java.util.Collection;
import java.util.Collections;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.LinkedList;
import java.util.Random;

import java.util.logging.Logger;
import java.util.logging.Level;

import java.util.concurrent.locks.*;
import java.util.concurrent.atomic.AtomicBoolean;

public class FingerprintScreener {
    static final Logger logger = 
	Logger.getLogger(FingerprintScreener.class.getName());

    private static final double LN2 = Math.log(2.);
    static final int DEFAULT_MAX_DEPTH = 10;

    
    static class SNode {
	int bit;
	SNode left, right;
	SNode parent;
	List<FPV> values = new ArrayList<FPV>();

	SNode () {}
	SNode (FPV v) { values.add(v); }
	SNode (List<FPV> values) {
	    this.values = values;
	}

	boolean isLeaf () { return left == null && right == null; }
    }

    static class FPV {
	Object key;
	int[] fp;
        int bits;

	FPV (Object key, int[] fp) {
	    this.key = key;
            for (int i = 0; i < fp.length; ++i) {
                bits += Integer.bitCount(fp[i]);
            }
	    this.fp = fp;
	}

	boolean get (int bit) { // in bit coordinate
	    return (fp[bit/32] & (1<<(bit%32))) != 0;
	}

	public String toString () {
	    StringBuilder sb = new StringBuilder ();
	    sb.append(key.toString()+"[");
	    sb.append(fp[0]+"="+Integer.toBinaryString(fp[0]));
	    for (int i = 1; i < fp.length; ++i) {
		sb.append(","+fp[i]+"="+Integer.toBinaryString(fp[i]));
	    }
	    sb.append("]");
	    return sb.toString();
	}
    }

    public static class Hit {
        Object key;
        double similarity;

        Hit (Object key, double similarity) {
            this.key = key;
            this.similarity = similarity;
        }

        public Object getKey () { return key; }
        public double getSimilarity () { return similarity; }
    }

    public static class ScreenStats {
	int hitCount;
	int screenCount;
	List<int[]> signatures = new ArrayList<int[]>();

	ScreenStats () {}
	ScreenStats (int hitCount, int screenCount) {
	    this.hitCount = hitCount;
	    this.screenCount = screenCount;
	}

	public int getHitCount () { return hitCount; }
	public int getScreenCount () { return screenCount; }
	public int[][] getSignatures () {
	    return signatures.toArray(new int[0][]);
	}
	public static int signaturePos (int sig) {
	    return sig >> 1;
	}
	public static int signatureBit (int sig) {
	    return sig & 1;
	}
    }

    public static class IndexStats {
	int nodeCount; // total node count = \sum_{k=0}^n 2^k
	int leafCount;
	int maxDepth;
	int minLeafSize, maxLeafSize;
	double avgLeafSize;
	int[] minsig, maxsig;
	int size;

	IndexStats () {}
	public int leafCount () { return leafCount; }
	public int maxDepth () { return maxDepth; }
	public int minLeafSize () { return minLeafSize; }
	public int maxLeafSize () { return maxLeafSize; }
	public double avgLeafSize () { return avgLeafSize; }
	public int size () { return size; }
	public int nodeCount () { return nodeCount; }
	public int[] minSignature () { return minsig; }
	public int[] maxSignature () { return maxsig; }

	public String toString () {
	    StringBuilder sb = new StringBuilder
		("##    NumNodes: "+nodeCount+"\n"
		 +"##    NumLeafs: "+leafCount+"\n"
		 +"##    MaxDepth: "+maxDepth+"\n"
		 +"## AvgLeafSize: "+avgLeafSize+"\n"
		 +"## MinLeafSize: "+minLeafSize+"\n"
		 +"## MaxLeafSize: "+maxLeafSize+"\n"
		 +"## ElementSize: "+size+"\n");

	    sb.append("##MinSignature:");
	    for (int i = 0; i < minsig.length; ++i) {
		sb.append(" ("+(minsig[i]>>1)+","+(minsig[i]&1)+")");
	    }
	    sb.append("\n##MaxSignature:");
	    for (int i = 0; i < maxsig.length; ++i) {
		sb.append(" ("+(maxsig[i]>>1)+","+(maxsig[i]&1)+")");
	    }
	    sb.append("\n");
	    return sb.toString();
	}
    }

    public enum ScreenType {
        IN, // default containtment
            EXACT, // exact match
            SUPER; // super match
    }

    private final ReentrantReadWriteLock rwl = new ReentrantReadWriteLock ();
    private final Lock rlock = rwl.readLock();
    private final Lock wlock = rwl.writeLock();
    private final AtomicBoolean isRebuilding = new AtomicBoolean (false);

    protected int dim;
    protected int maxdepth;
    protected SNode root = null;

    volatile protected int[] _prof;
    volatile protected List<FPV> values = new ArrayList<FPV>();

    public FingerprintScreener (int dim) {
	this (dim, DEFAULT_MAX_DEPTH);
    }

    public FingerprintScreener (int dim, int maxdepth) {
	if (dim < 32) {
	    throw new IllegalArgumentException
		("Bad dimension (< 32) specified: "+dim);
	}
	else if ((dim % 32) != 0) {
	    throw new IllegalArgumentException
		("Dimension "+dim+" is not a multiple of 32!");
	}
	this.dim = dim;
	this.maxdepth = maxdepth;

	_prof = new int[dim];
    }

    public void add (Object key, int[] fp) {
	if (isRebuilding.get()) {
	    throw new IllegalStateException
		("Index is currently being rebuilt");
	}

	if (fp.length*32 != dim) {
	    throw new IllegalArgumentException
		("Entry "+key+" has invalid dimension: "+(fp.length*32));
	}
	
	wlock.lock();
	try {
	    addFast (key, fp);
	}
	finally {
	    wlock.unlock();
	}
    }

    /**
     * Use this only if you know for sure only one thread is accessing
     * this instance and that the fingerprint dimension is correct!
     */
    public void addFast (Object key, int[] fp) {
	FPV val = new FPV (key, fp);
	if (root == null) {
	    root = new SNode (val);
	}
	else {
	    insert (val);
	}
	values.add(val);
    }

    public int size () { 
	rlock.lock();
	try {
	    return values.size(); 
	}
	finally {
	    rlock.unlock();
	}
    }

    public int getDim () { return dim; }

    /**
     * Rebuilding the index with a different maxdepth parameter.  When this
     * method is called, it blocks all other methods (e.g., search, add).
     */
    public synchronized void rebuild (int maxdepth) {
	if (rwl.hasQueuedThreads()) {
	    throw new IllegalStateException 
		("Index is currently being used; "
                 +"can't rebuild at the moment!");
	}

	isRebuilding.set(true);
        try {
            this.maxdepth = maxdepth;
            
            Iterator<FPV> iter = values.iterator();
            root = new SNode (iter.next());
            while (iter.hasNext()) {
		insert (iter.next());
	    }
	}
        finally {
            isRebuilding.set(false);
        }
    }

    static String toString (SNode node) {
	StringBuilder sb = new StringBuilder ();
	depthFirst (sb, 0, node, "**");
	return sb.toString();
    }

    static void depthFirst (StringBuilder sb, int depth, 
			    SNode n, String prefix) {
	if (n == null) {
	    return;
	}
	
	for (int i = 0; i <= depth; ++i) {
	    sb.append("  ");
	}

	if (prefix != null) {
	    sb.append(prefix);
	}
	if (n.isLeaf()) {
	    sb.append(" d="+depth);
	    for (FPV v : n.values) {
		sb.append(" "+v);
	    }
	    sb.append("\n");
	}
	else {
	    sb.append(" d="+depth+" bit="+n.bit+"\n");
	}
	depthFirst (sb, depth+1, n.left, "<");
	depthFirst (sb, depth+1, n.right,">");
    }

    private void insert (FPV x) {
	LinkedList<SNode> stack = new LinkedList<SNode>();
	stack.push(root);

	int depth = 0;
	while (!stack.isEmpty()) {
	    SNode v = stack.pop();
	    if (v.isLeaf()) {
		int k = dim;

		if (maxdepth <= 0 || depth < maxdepth) {
		    for (FPV w : v.values) {
			for (int i = 0; i < dim; ++i) {
			    if (w.get(i) != x.get(i)) {
				++_prof[i];
			    }
			}
		    }
		    int max = 0;
		    for (int i = 0; i < _prof.length; ++i) {
			if (_prof[i] > max) {
			    max = _prof[i];
			    k = i;
			}
			_prof[i] = 0;
		    }
		}

		if (k == dim) {
		    v.values.add(x);
		}
		else {
		    v.bit = k; // promote this leaf into internal node
		    if (x.get(k)) {
			v.right = new SNode (x);
			v.left = new SNode (v.values);
		    }
		    else {
			v.right = new SNode (v.values);
			v.left = new SNode (x);
		    }
		    v.right.parent = v;
		    v.left.parent = v;
		    v.values = null;
		}
	    }
	    else {
		SNode child = x.get(v.bit) ? v.right : v.left;
		stack.push(child);
		++depth;
	    }
	}
    }

    SNode createNode (List<FPV> subset, int depth) {
	if (subset == null || subset.isEmpty()) {
	    return null;
	}

	int[] freq = new int[dim];
	for (FPV v : subset) {
	    for (int i = 0; i < dim; ++i) {
		if (v.get(i)) {
		    ++freq[i];
		}
	    }
	}

	// find an index that split the set into as evenly in two
	// halves as possible
	int min = Integer.MAX_VALUE, split = -1, half = subset.size()/2;
	for (int i = 0; i < dim; ++i) {
	    int d = Math.abs(freq[i]-half);
	    if (d < min) {
		min = d;
		split = i;
	    }
	}

	SNode node = new SNode (subset);
	node.bit = split;
	if (min == half || split < 0 || depth >= maxdepth) {
	    // leaf node
	}
	else {
	    List<FPV> lsub = new ArrayList<FPV>();
	    List<FPV> rsub = new ArrayList<FPV>();
	    for (FPV v : subset) {
		if (v.get(split)) {
		    rsub.add(v);
		}
		else {
		    lsub.add(v);
		}
	    }
	    
	    node.left = createNode (lsub, depth+1);
	    if (node.left != null) {
		node.left.parent = node;
	    }

	    node.right = createNode (rsub, depth+1);
	    if (node.right != null) {
		node.right.parent = node;
	    }
	}

	return node;
    }

    public ScreenStats screen (int[] fp) {
	return screen (fp, ScreenType.IN, null);
    }

    public ScreenStats screen (int[] fp, ScreenType type) {
	return screen (fp, type, null);
    }

    int[] getSignature (SNode n) {
	List<Integer> sig = new ArrayList<Integer>();
	for (SNode p = n.parent; p != null; p = p.parent) {
	    sig.add((p.bit << 1) | (p.left == n ? 0 : 1));
	    if (p.left != n && p.right != n) {
		System.err.println("FATAL ERROR!");
		System.exit(1);
	    }
	    n = p;
	}
	int[] ps = new int[sig.size()];
	for (int i = 0; i < ps.length; ++i) { ps[i] = sig.get(i); }
	return ps;
    }

    public ScreenStats screen (int[] fp, ScreenType type,
                               SearchCallback<Hit> callback) {
	if (isRebuilding.get()) {
	    throw new IllegalStateException
		("Index is currently being rebuilt!");
	}

	if (fp.length*32 != dim) {
	    throw new IllegalArgumentException
		("Query has invalid dimension: "+(fp.length*32));
	}
	if (root == null) {
	    throw new IllegalArgumentException
		("Signature tree hasn't been constructed yet!");
	}

	rlock.lock();
	try {
            switch (type) {
            case SUPER:
                return screenSUPER (fp, callback);

            case EXACT:
                return screenEXACT (fp, callback);

            case IN:
            default:
                return screenIN (fp, callback);
            }
	}
	finally {
	    rlock.unlock();
	}
    }

    ScreenStats screenSUPER (int[] fp, SearchCallback<Hit> callback) {
        LinkedList<SNode> stack = new LinkedList<SNode>();
        stack.push(root);
	
        ScreenStats stats = new ScreenStats ();
        while (!stack.isEmpty()) {
            SNode n = stack.pop();
            if (n.isLeaf()) {
                int[] sig = getSignature (n);
                stats.signatures.add(sig);
                
                if (!screenSUPER (stats, n.values, fp, callback)) {
                    break;
                }
            }
            else {
                if ((fp[n.bit/32] & (1<<(n.bit%32))) != 0) {
                    stack.push(n.right);
                }
                stack.push(n.left);
            }
        }
        return stats;
    }

    ScreenStats screenEXACT (int[] fp, SearchCallback<Hit> callback) {
        LinkedList<SNode> stack = new LinkedList<SNode>();
        stack.push(root);
	
        ScreenStats stats = new ScreenStats ();
        while (!stack.isEmpty()) {
            SNode n = stack.pop();
            if (n.isLeaf()) {
                int[] sig = getSignature (n);
                stats.signatures.add(sig);
                
                if (!screenEXACT (stats, n.values, fp, callback)) {
                    break;
                }
            }
            else {
                stack.push(n.right);
                if ((fp[n.bit/32] & (1<<(n.bit%32))) == 0) {
                    stack.push(n.left);
                }
            }
        }
        return stats;
    }


    ScreenStats screenIN (int[] fp, SearchCallback<Hit> callback) {
        LinkedList<SNode> stack = new LinkedList<SNode>();
        stack.push(root);
	
        ScreenStats stats = new ScreenStats ();
        while (!stack.isEmpty()) {
            SNode n = stack.pop();
            if (n.isLeaf()) {
                int[] sig = getSignature (n);
                stats.signatures.add(sig);
                
                if (!screenIN (stats, n.values, fp, callback)) {
                    break;
                }
            }
            else {
                stack.push(n.right);
                if ((fp[n.bit/32] & (1<<(n.bit%32))) == 0) {
                    stack.push(n.left);
                }
            }
        }
        return stats;
    }
    

    boolean screenIN (ScreenStats stats, Collection<FPV> values, 
                      int[] fp, SearchCallback<Hit> callback) {
	for (FPV v : values) {
	    int i = 0;
	    for (; i < fp.length; ++i) {
		if ((v.fp[i] & fp[i]) != fp[i]) {
		    break;
		}
	    }
	    ++stats.screenCount;
	    
	    boolean matched = i == fp.length;
	    if (matched) {
		if (callback != null) {
                    Hit h = new Hit (v.key, calcTanimoto (fp, v.fp));
		    if (!callback.matched(h)) {
			return false;
		    }
		}
		++stats.hitCount;
	    }
	}

	return true;
    }

    boolean screenSUPER (ScreenStats stats, Collection<FPV> values, 
                         int[] fp, SearchCallback<Hit> callback) {
	for (FPV v : values) {
	    int i = 0;
	    for (; i < fp.length; ++i) {
		if ((v.fp[i] & fp[i]) != v.fp[i]) {
		    break;
		}
	    }
	    ++stats.screenCount;
	    
	    boolean matched = i == fp.length;
	    if (matched) {
		if (callback != null) {
                    Hit h = new Hit (v.key, calcTanimoto (fp, v.fp));
		    if (!callback.matched(h)) {
			return false;
		    }
		}
		++stats.hitCount;
	    }
	}

	return true;
    }

    boolean screenEXACT (ScreenStats stats, Collection<FPV> values, 
                         int[] fp, SearchCallback<Hit> callback) {
	for (FPV v : values) {
	    int i = 0;
	    for (; i < fp.length; ++i) {
		if (v.fp[i] != fp[i]) {
		    break;
		}
	    }
	    ++stats.screenCount;
	    
	    boolean matched = i == fp.length;
	    if (matched) {
		if (callback != null) {
                    Hit h = new Hit (v.key, calcTanimoto (fp, v.fp));
		    if (!callback.matched(h)) {
			return false;
		    }
		}
		++stats.hitCount;
	    }
	}

	return true;
    }


    /**
     * threshold - specify the minimum similarity value
     */
    public ScreenStats screen (int[] fp, double threshold) {
        return screen (fp, threshold, null);
    }

    public ScreenStats screen (int[] fp, double threshold, 
                               SearchCallback<Hit> callback) {
	if (fp.length*32 != dim) {
	    throw new IllegalArgumentException
		("Query has invalid dimension: "+(fp.length*32));
	}

        if (threshold < 0.001) {
            throw new IllegalArgumentException
                ("Similarity threshold is too small: "+threshold);
        }

        int bits = 0;
        for (int i = 0; i < fp.length; ++i) {
            bits += Integer.bitCount(fp[i]);
        }

        ScreenStats stats = new ScreenStats ();
        for (FPV v : values) {
            if (v.bits > 0) {
                // do a quick bound checking...
                double t = (double)Math.min(bits, v.bits)
                    / Math.max(bits, v.bits);
                if (threshold <= t) {
                    t = calcTanimoto (fp, v.fp);
                    if (t >= threshold) {
                        if (callback != null) {
                            Hit h = new Hit (v.key, t);
                            if (!callback.matched(h)) {
                                break;
                            }
                        }
                        ++stats.hitCount;
                    }
                }
            }
            ++stats.screenCount;
        }

        return stats;
    }

    static double calcTanimoto (int[] fp1, int[] fp2) {
        int a = 0, b = 0;
        for (int i = 0; i < fp1.length; ++i) {
            a += Integer.bitCount(fp1[i] & fp2[i]);
            b += Integer.bitCount(fp1[i] | fp2[i]);
        }
        return b == 0 ? 0. : (double)a/b;
    }


    /*
     * perform linear scan
     */
    public int linear (int[] fp) { 
	return linear (fp, null);
    }

    public int linear (int[] fp, SearchCallback callback) {
	int i, count = 0;
	for (FPV v : values) {
	    i = 0;
	    for (; i < fp.length; ++i) {
		if ((v.fp[i] & fp[i]) != v.fp[i]) {
		    break;
		}
	    }

	    boolean matched = i == fp.length;
	    if (matched) {
		if (callback != null) {
		    callback.matched(v.key);
		}
		++count;
	    }
	}
	return count;
    }

    public void dump (PrintStream ps) {
	ps.println(toString (root));
    }

    public IndexStats getIndexStats () {
	if (isRebuilding.get()) {
	    throw new IllegalStateException
		("Index is currently being rebuilt!");
	}

	rlock.lock();
	try {
	    return _getIndexStats ();
	}
	finally {
	    rlock.unlock();
	}
    }

    IndexStats _getIndexStats () {
	if (root == null) {
	    return null;
	}

	LinkedList<SNode> stack = new LinkedList<SNode>();
	stack.push(root);

	IndexStats stats = new IndexStats ();
	stats.minLeafSize = Integer.MAX_VALUE;

	while (!stack.isEmpty()) {
	    SNode n = stack.pop();
	    if (n.isLeaf()) {
		int d = stack.size();
		if (d > stats.maxDepth) {
		    stats.maxDepth = d;
		}

		int size = n.values.size();
		if (size < stats.minLeafSize) {
		    stats.minLeafSize = size;
		    stats.minsig = getSignature (n);
		}
		if (size > stats.maxLeafSize) {
		    stats.maxLeafSize = size;
		    stats.maxsig = getSignature (n);
		}
		stats.avgLeafSize += size;

		++stats.leafCount;
	    }
	    else {
		stack.push(n.right);
		stack.push(n.left);
	    }
	    ++stats.nodeCount;
	}
	stats.avgLeafSize /= stats.leafCount;
	stats.size = values.size();

	return stats;
    }


    static void testRandomScreen (int dim, int size, int test) {
	java.util.Random rand = new java.util.Random (1l);
	FingerprintScreener screener = new FingerprintScreener (dim, 12);

	logger.info("generating "+size+" random "
		    +dim+"-bit fingerprints...");

	for (int i = 0; i < size; ++i) {
	    int[] fp = randomFp (rand, dim, rand.nextGaussian());
	    screener.add("key"+i, fp);
	}
	FingerprintScreener.IndexStats stats = screener.getIndexStats();
	logger.info("** Signature tree stats\n"+stats);

	logger.info("testing screening performance...");
	for (int i = 0; i < test; ++i) {
	    int[] fp = randomFp (rand, dim);
	    System.out.print("Query: ["+fp[0]);
	    for (int j = 1; j< fp.length; ++j) {
		System.out.print(","+fp[j]);
	    }
	    System.out.println("]");
	    double den = 0.;
	    for (int j = 0; j < fp.length; ++j) {
		den += Integer.bitCount(fp[j]);
	    }
	    den /= dim;
	    System.out.println
		("Query density: "+String.format("%1$.3f", den));

            final List<Hit> hits = new ArrayList<Hit>();
	    long start = System.currentTimeMillis();
	    ScreenStats ss = screener.screen
                (fp, 0.6, new SearchCallback<Hit> () {
                    public boolean matched (Hit h) {
                        hits.add(h);
                        return true;
                    }
                });
	    double time = 1e-3*(System.currentTimeMillis()-start);
	    System.out.println("Signature screen took "
			       +String.format("%1$.3fs",time)
			       +" spanning "+ss.getSignatures().length
			       +" bucket(s); "+ss.getHitCount()+"/"
			       +ss.getScreenCount()+" found!");
	    start = System.currentTimeMillis();
	    int cnt2 = screener.linear(fp);
	    time = 1e-3*(System.currentTimeMillis()-start);
	    System.out.println("Linear screen took "
			       +String.format("%1$.3fs",time)
			       +"; "+cnt2+" found!");
            /*
	    if (ss.getHitCount() != cnt2) {
		System.out.println("** FATAL: bug found for this query; "
				   +"please report this error!");
		System.exit(1);
	    }
            */
            if (!hits.isEmpty()) {
                double s = 0.;
                for (Hit h : hits) {
                    s += h.getSimilarity();
                }
                s /= hits.size();
                System.out.println("Avg similarity: "
                                   +String.format("%1$.3f", s));
            }
	    System.out.println();
	}	
    }

    static int[] randomFp (java.util.Random rand, int dim) {
	return randomFp (rand, dim, rand.nextDouble());
    }

    static int[] randomFp (java.util.Random rand, int dim, double density) {
	int[] fp = new int[dim/32];
	if (density < 0) density *= -1.;
	int nb = (int)(density*dim + .5);
	for (int i = 0; i < nb; ++i) {
	    int b = rand.nextInt(dim); // uniformly turn on the bits
	    fp[b/32] |= 1<<(b%32);
	}
	return fp;
    }

    static void test () {
	FingerprintScreener fps = new FingerprintScreener (32);
	Random rand = new Random ();
	for (int i = 0; i < 20; ++i) {
	    int[] fp = new int[1];
	    fp[0] = rand.nextInt(20);
	    fps.add("key"+i, fp);
	}
	System.out.println("** Signature Tree **\n");
	fps.dump(System.out);
	IndexStats stats = fps.getIndexStats();
	logger.info("** Signature tree stats\n"+stats);	
    }

    public static void main (String[] argv) throws Exception {
	testRandomScreen (1024, 1000000, 100);
	//test ();
    }
}
