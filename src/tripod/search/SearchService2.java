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

import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicLong;
import java.util.logging.Logger;
import java.util.logging.Level;

import chemaxon.sss.SearchConstants;
import chemaxon.sss.search.*;
import chemaxon.struc.Molecule;
import chemaxon.struc.MolAtom;
import chemaxon.util.MolHandler;


public class SearchService2 implements MoleculeService {
    private static final Logger logger = 
	Logger.getLogger(SearchService2.class.getName());

    protected static final int DEFAULT_PARTITION_COUNT = 2;
    protected static final int DEFAULT_MAXDEPTH = 10;
    protected static final int DEFAULT_FP_SIZE = 16; // 16*32 = 512 bits
    protected static final int DEFAULT_FP_DARKNESS = 2; // 2 bits for each pattern
    protected static final int DEFAULT_FP_DEPTH = 6; // recursion depth


    static public class MolEntry implements Comparable<MolEntry> {
        Object key;
        Molecule mol;
        int[] fingerprint;
        int[][] mappings;
        Double rank;
        double similarity;

        MolEntry (Molecule mol) {
            this.mol = mol;
        }

        MolEntry (Object key, Molecule mol) {
            this.key = key;
            this.mol = mol;
        }

        MolEntry (Object key, Molecule mol, int[] fingerprint) {
            this.key = key;
            this.mol = mol;
            this.fingerprint = fingerprint;
        }

        public Molecule getMol () { return mol; }
        public Object getKey () { return key; }

        public int[] getFingerprint () { return fingerprint; }
        public void setFingerprint (int[] fingerprint) { 
            this.fingerprint = fingerprint; 
        }
        
        public double getSimilarity () { return similarity; }
        public void setSimilarity (double similarity) { 
            this.similarity = similarity; 
        }

        public Double getRank () { return rank; }
        public void setRank (Double rank) { this.rank = rank; }

        public void setAtomMappings (int[][] mappings) {
            this.mappings = mappings;
        }
        public int[][] getAtomMappings () { return mappings; }

        public int compareTo (MolEntry me) {
            if (rank != null) {
                return rank.compareTo(me.rank);
            }
            if (rank != null || me.rank != null)
                return 1;
            return 0;
        }
    }

    static class Query {
        MolEntry query;
        SearchParams params;

        public Query (MolEntry query) {
            this (query, SearchParams.substructure());
        }

        public Query (MolEntry query, SearchParams params) {
            this.query = query;
            this.params = params;
        }

        public MolEntry getEntry () { return query; }
        public SearchParams getParams () { return params; }
    }

    class SearchWorker implements Callable<Integer>, 
                                  SearchCallback<FingerprintScreener.Hit> {
        Query query;
        FingerprintScreener screener;
        SearchCallback<MolEntry> callback;
        final CountDownLatch signal;
        int matches;

        MolSearch msearch = new MolSearch ();
        Molecule qmol;

        /*
         * if callback is null, then we simply do counting
         */
        SearchWorker (Query query, FingerprintScreener screener,
                      SearchCallback<MolEntry> callback,
                      CountDownLatch signal) {
            this.query = query;
            this.signal = signal;
            this.callback = callback;
            this.screener = screener;

            this.qmol = query.getEntry().getMol().cloneMolecule();
        }

        /*
         * Callable
         */
        public Integer call () {
            FingerprintScreener.ScreenStats stats = null;
            try {
                int[] fp = query.getEntry().getFingerprint();

                switch (query.getParams().getType()) {
                case Superstructure:
                    stats = screener.screen
                        (fp, FingerprintScreener.ScreenType.SUPER, 
                         callback != null ? this : null);
                    break;

                case Substructure:
                    stats = screener.screen
                        (fp, FingerprintScreener.ScreenType.IN, 
                         callback != null ? this : null);
                    break;

                case Similarity:
                    stats = screener.screen
                        (fp, query.getParams().getSimilarity(), 
                         callback != null ? this : null);
                    break;

                case Exact:
                    stats = screener.screen
                        (fp, FingerprintScreener.ScreenType.EXACT, 
                         callback != null ? this : null);
                }
            }
            finally {
                if (signal != null)
                    signal.countDown();
            }

            return callback != null ? matches : stats.getHitCount();
        }
        
        /*
         * SearchCallback
         */
        public boolean matched (FingerprintScreener.Hit hit) {
            Molecule mol = getMol (hit.getKey());
            if (mol == null) {
                return true; // continue
            }

            Molecule target = mol.cloneMolecule();
            target.aromatize();

            MolEntry entry = null;
            switch (query.getParams().getType()) {
            case Substructure: 
                entry = doSubstructure (hit, target); 
                break;

            case Superstructure:
                entry = doSuperstructure (hit, target);
                break;

            case Similarity:
                entry = doSimilarity (hit, target);
                break;

            case Exact:
                entry = doExact (hit, target);
            }

            boolean truncated = false;
            if (entry != null && callback != null) {
                truncated = !callback.matched(entry);
            }

	    if (matches < query.getParams().getMatchLimit()) {
	    }
	    else {
		logger.warning("** Query "
                               + qmol.toFormat("cxsmarts") +
			       " exceeds match limit "+matches+"/"
			       + query.getParams().getMatchLimit());
		truncated = true;
	    }
            
            return !Thread.currentThread().isInterrupted() && !truncated;
        }

        MolEntry doSubstructure (FingerprintScreener.Hit hit, 
                                 Molecule target) {
            // this assumes that the query mol has already been 
            // properly aromatize (via the fingerprint generation)
            msearch.setSearchType(SearchConstants.SUBSTRUCTURE);
	    msearch.setQuery(qmol);
	    msearch.setTarget(target);

            return doSearch (hit, target);
        }

        MolEntry doSuperstructure (FingerprintScreener.Hit hit, 
                                   Molecule target) {
            msearch.setSearchType(SearchConstants.SUPERSTRUCTURE);
            msearch.setQuery(qmol);
            msearch.setTarget(target);

            return doSearch (hit, target);
        }

        MolEntry doSimilarity (FingerprintScreener.Hit hit, Molecule target) {
            target.dearomatize();
            MolEntry entry = new MolEntry (hit.getKey(), target);
            entry.setSimilarity(hit.getSimilarity());
            ++matches;

            return entry;
        }

        MolEntry doExact (FingerprintScreener.Hit hit, Molecule target) {
            msearch.setSearchType(SearchConstants.EXACT);
            msearch.setQuery(qmol);
            msearch.setTarget(target);

            return doSearch (hit, target);
        }

        MolEntry doSearch (FingerprintScreener.Hit hit, Molecule target) {

            MolEntry entry = null;
	    try {
		int[][] hits = msearch.findAll();
		if (hits != null) {
		    ++matches;
		    
		    target.dearomatize();
                    entry = new MolEntry (hit.getKey(), target);
                    entry.setAtomMappings(hits);
                    entry.setSimilarity(hit.getSimilarity());
		}
	    }
	    catch (SearchException ex) {
		logger.log(Level.SEVERE, 
			   "Search fails for query "
			   +qmol.toFormat("cxsmarts"), ex);
	    }

            return entry;
        }
    } // SearchWorker

    public static int[] generateFingerprint (Molecule mol, int dim) {
        return generateFingerprint (mol, dim, true);
    }

    public static int[] generateFingerprint 
        (Molecule mol, int dim, boolean aromatize) {
	MolHandler mh = new MolHandler (mol);
        if (aromatize) { // aromatize
            mh.aromatize();
        }
	return mh.generateFingerprintInInts
	    (dim, DEFAULT_FP_DARKNESS, DEFAULT_FP_DEPTH);
    }

    public static Molecule createMol (String str) {
	try {
	    MolHandler mh = new MolHandler (str, true);
	    Molecule mol = mh.getMolecule();
            //mol.hydrogenize(false);
            return mol;
	}
	catch (Exception ex) {
	    throw new IllegalArgumentException ("Invalid molecule");
	}
    }

    // fingerprint dimension in int's
    protected int dim;
    // current size of all indexes
    protected AtomicLong size = new AtomicLong(); 
    protected long capacity; // maximum index size; 0 => no limits
    protected volatile FingerprintScreener[] screeners;
    protected ExecutorService threadPool = Executors.newCachedThreadPool();

    public SearchService2 () {
        this (DEFAULT_PARTITION_COUNT, DEFAULT_FP_SIZE);
    }

    public SearchService2 (int partitions) {
        this (partitions, DEFAULT_FP_SIZE);
    }

    public SearchService2 (int partitions, int dim) {
        this.dim = dim;
        setPartitionCount (partitions);
    }

    public void add (Object key, int[] fingerprint) {
        int index = getPartition (key);
        screeners[index].add(key, fingerprint);
        size.incrementAndGet();
    }

    public FingerprintScreener[] getIndexes () { return screeners; }

    public void setThreadPool (ExecutorService threadPool) {
        try {
            this.threadPool.shutdown();
        }
        finally {
            this.threadPool = threadPool;
        }
    }
    public ExecutorService getThreadPool () { return threadPool; }

    public void setPartitionCount (int partitions) {
        if (partitions <= 0) {
            throw new IllegalArgumentException ("Partition size must be > 0");
        }

        if (screeners != null && partitions == screeners.length) {
            return; // same partition
        }

        screeners = new FingerprintScreener[partitions];
        for (int i = 0; i < partitions; ++i) {
            screeners[i] = new FingerprintScreener (dim*32);
        }
        size.set(0);
    }
    public int getPartitionCount () { return screeners.length; }

    public void clear () {
        screeners = new FingerprintScreener[screeners.length];
        for (int i = 0; i < screeners.length; ++i) {
            screeners[i] = new FingerprintScreener (dim*32);
        }
        size.set(0);
    }

    public void setCapacity (long capacity) { 
        this.capacity = capacity;
    }
    public long getCapacity () { return capacity; }
    public long size () { return size.get(); }

    public int count (String mol) {
        return search (createMol (mol), SearchParams.substructure(), null);
    }
    public int count (String mol, SearchParams params) {
        return search (createMol (mol), params, null);
    }
    public int count (Molecule mol) {
        return search (mol, SearchParams.substructure(), null);
    }
    public int count (Molecule mol, SearchParams params) {
        return search (mol, params, null);
    }
    protected int count (Query query) {
        return search (query, null);
    }

    public int search (String mol, SearchCallback<MolEntry> callback) {
        return search (createMol (mol), SearchParams.substructure(), callback);
    }
    public int search (String mol, SearchParams params, 
                       SearchCallback<MolEntry> callback) {
        return search (createMol (mol), params, callback);
    }
    public int search (Molecule mol, SearchCallback<MolEntry> callback) {
        return search (mol, SearchParams.substructure(), callback);
    }

    public int search (Molecule mol, SearchParams params, 
                       SearchCallback<MolEntry> callback) {
        MolEntry entry = new MolEntry (mol);
        entry.setFingerprint(generateFingerprint
                             (mol, getDim (), params.getAromatize()));

        double den = calcDensity (entry.getFingerprint());
        if (den < params.getMinDensity()) {
            logger.warning("Query doesn't have sufficient density!");
            return -1;
        }

        return search (new Query (entry, params), callback);
    }


    protected int search (Query query, SearchCallback<MolEntry> callback) {
        if (screeners.length == 1) {
            return new SearchWorker 
                (query, screeners[0], callback, null).call();
        }

        CountDownLatch signal = new CountDownLatch (screeners.length);

        List<Future<Integer>> tasks = new ArrayList<Future<Integer>>();
        for (FingerprintScreener s : screeners) {
            Future<Integer> f = threadPool.submit
                (new SearchWorker (query, s, callback, signal));
            tasks.add(f);
        }

        int total = 0;
        try {
            long timeout = query.getParams().getTimeout();
            if (timeout > 0) {
                if (!signal.await(timeout, TimeUnit.MILLISECONDS)) {
                    for (Future<Integer> f : tasks) {
                        try {
                            f.cancel(true);
                        }
                        catch (Exception ex) {
                            ex.printStackTrace();
                        }
                    }

                    logger.warning
                        ("Search time expired: "+timeout+"ms");
                }
            }
            else {
                signal.await();
            }
        }
        catch (InterruptedException ex) {
            logger.log(Level.SEVERE, 
                       "Thread "+Thread.currentThread()+" interrupted!", ex);
        }
        finally {
            for (Future<Integer> f : tasks) {
                try {
                    Integer c = f.get();
                    if (c != null) {
                        total += c;
                    }
                }
                catch (Exception ex) {
                    ex.printStackTrace();
                }
            }
            //logger.info(total+" matches!");
        }

        return total;
    }

    static double calcDensity (int[] fp) {
        double c = 0.;
        for (int i = 0; i < fp.length; ++i) {
            c += Integer.bitCount(fp[i]);
        }
        c /= 32*fp.length;
        return c;
    }

    public void shutdown () {
        try {
            threadPool.shutdownNow();
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    /**
     * A place holder method to allow subclass to initialize
     * this search service.
     */
    public void init (Object ctx) {
    }

    public int getDim () { return dim; }
    public void setDim (int dim) {
        if (size.get() > 0) {
            throw new IllegalStateException
                ("Can't set dim once index is built");
        }
        this.dim = dim;
    }

    /**
     * Sub class should provide a mechanism to retrieve a Molecule
     * instance given a key object that was used during the construction
     * of the indexes.  This method must be thread safe as it's called
     * from multiple threads!
     */
    public Molecule getMol (Object key) {
	if (key instanceof Molecule) {
	    return (Molecule)key;
	}
	return null;
    }

    /**
     * Subclass can override this to provide a custom means of cluster
     * the keys. This method must be thread-safe!
     */
    protected int getPartition (Object key) {
        long h = key.hashCode() & 0xffffffffl;
        return (int)(h % screeners.length);
    }
}
