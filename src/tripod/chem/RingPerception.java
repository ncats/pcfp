package tripod.chem;

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
import java.util.logging.Logger;
import java.util.logging.Level;

import java.io.InputStream;
import java.io.FileInputStream;
import java.io.File;

import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;
import chemaxon.formats.MolImporter;


/**
 * A simple ring perception implementation based on 
 * J.D. Horton, A polynomial-time algorithm to find the shortest cycle
 *   basis of a graph. SIAM J. Comput., 16(2):358-366, 1987.
 */
public class RingPerception {
    private static final Logger logger = 
	Logger.getLogger(RingPerception.class.getName());

    /**
     * default max ring size
     */
    static final int MAX_RING_SIZE = 20;

    /**
     * default max number of rings to generate
     */
    static final int MAX_NUM_RINGS = 128;

    static class Path {
        BitSet aset = new BitSet ();
        BitSet bset = new BitSet ();
        int[] atoms, bonds;

        MolAtom[] atomArray;
        MolBond[] bondArray;
        Molecule mol;

        Path (Collection<MolBond> path) {
            add (path.toArray(new MolBond[0]));
        }

        Path (MolBond... path) {
            add (path);
        }

        Path (Path... paths) {
            List<MolBond> bonds = new ArrayList<MolBond>();
            for (Path p : paths) {
                for (MolBond b : p.getBondArray()) 
                    bonds.add(b);
            }

            add (bonds.toArray(new MolBond[0]));
        }

        protected void add (Path path) {
            add (path.getBondArray());
        }

        public void add (MolBond...path) {
            int i = 0;
            for (MolBond b : path) {
                if (mol == null) {
                    mol = (Molecule)b.getParent();
                }
                else if (mol != b.getParent()) {
                    throw new IllegalArgumentException 
                        ("Input path contains bond from "
                         +"different parent molecules!");
                }

                int idx = mol.indexOf(b.getAtom1());
                aset.set(idx);

                idx = mol.indexOf(b.getAtom2());
                aset.set(idx);

                idx = mol.indexOf(b);
                bset.set(idx);
            }

            i = 0;
            atoms = new int[aset.cardinality()];
            atomArray = new MolAtom[atoms.length];
            for (int j = aset.nextSetBit(0); j >= 0; 
                 j = aset.nextSetBit(j+1)) {
                atoms[i] = j;
                atomArray[i] = mol.getAtom(j);
                ++i;
            }

            i = 0;
            bonds = new int[bset.cardinality()];
            bondArray = new MolBond[bonds.length];
            for (int j = bset.nextSetBit(0); j >= 0; 
                 j = bset.nextSetBit(j+1)) {
                bonds[i] = j;
                bondArray[i] = mol.getBond(j);
                ++i;
            }
        }

        public int hashCode () { return aset.hashCode(); }
        public boolean equals (Object obj) {
            if (obj instanceof Path) {
                return aset.equals(((Path)obj).aset);
            }
            return false;
        }

        public boolean intersects (Path path) {
            return aset.intersects(path.aset);
        }

        protected boolean overlaps (Path path, int atom) {
            if (aset.intersects(path.aset)) {
                BitSet bs = (BitSet)aset.clone();
                bs.and(path.aset);
                return bs.cardinality() == 1 && bs.nextSetBit(0) == atom;
            }
            return false;
        }

        public boolean hasAtom (int atom) { return aset.get(atom); }
        public boolean hasBond (int bond) { return bset.get(bond); }
        public BitSet getAtomSet () { return aset; }
        public BitSet getBondSet () { return bset; }
        public int[] getAtoms () { return atoms; }
        public int[] getBonds () { return bonds; }
        public int size () { return atoms.length; }
        public int getAtom (int i) { return atoms[i]; }
        public MolBond[] getBondArray () { return bondArray; }
        public MolAtom[] getAtomArray () { return atomArray; }
        public MolBond getBond (int a1, int a2) {
            if (aset.get(a1) && aset.get(a2)) {
                for (int i = 0; i < bonds.length; ++i) {
                    MolBond b = mol.getBond(bonds[i]);
                    int k = mol.indexOf(b.getAtom1());
                    int j = mol.indexOf(b.getAtom2());
                    if ((k == a1 && j == a2) || (k == a2 && j == a1))
                        return b;
                }
            }
            return null;
        }
        public boolean isConnected (int a1, int a2) {
            return getBond (a1, a2) != null;
        }

        public String toString () {
            return "{size="+size()+",atoms="+aset+",bonds="+bset+"}";
        }
    }

    static public class Ring extends Path {
        Ring (Path... paths) {
            super (paths);

            int[] degree = new int[mol.getAtomCount()];
            for (MolBond b : getBondArray ()) {
                int a1 = mol.indexOf(b.getAtom1());
                int a2 = mol.indexOf(b.getAtom2());
                ++degree[a1];
                ++degree[a2];
            }

            for (int i = 0; i < atoms.length; ++i) 
                if (degree[atoms[i]] != 2) {
                    throw new IllegalArgumentException
                        ("Specified paths do not form a cycle!");
                }
        }
    }

    private Molecule mol;
    private int[][] cost;
    private int[][] btab;
    private MolAtom[] atoms;
    private Set<Ring> rings = new HashSet<Ring>();
    private int maxsize;
    private int maxrings = MAX_NUM_RINGS;
    
    public RingPerception () {
        this (MAX_RING_SIZE);
    }

    public RingPerception (int maxsize) {
        this.maxsize = maxsize;
    }

    public RingPerception (Molecule mol) {
        this (mol, MAX_RING_SIZE);
    }

    public RingPerception (Molecule mol, int maxsize) {
        this (maxsize);
        perceive (mol);
    }

    public int getMaxSize () { return maxsize; }
    public void setMaxSize (int maxsize) { this.maxsize = maxsize; }

    public int getMaxRings () { return maxrings; }
    public void setMaxRings (int maxrings) { this.maxrings = maxrings; }

    public RingPerception perceive (Molecule mol) {
        this.mol = mol;
        atoms = mol.getAtomArray();

        int acount = atoms.length;
        cost = new int[acount][acount];
        btab = mol.getBtab();

        // init the path cost...
        for (int i = 0; i < acount; ++i) {
            cost[i][i] = 0;
            for (int j = i+1; j < acount; ++j) {
                cost[i][j] = cost[j][i] = btab[i][j] < 0 ? acount : 1;
            }
        }
        
        // 1. now perform floyd-warshall's all pairs shortest path
        for (int k = 0; k < acount; ++k) 
            for (int i = 0; i < acount; ++i) 
                for (int j = 0; j < acount; ++j) 
                    cost[i][j] = Math.min
                        (cost[i][j], cost[i][k]+cost[k][j]);

        rings.clear();

        // 2. find all cycles of at most <= maxsize
        MolBond[] bonds = mol.getBondArray();
        for (int i = 0; i < atoms.length; ++i)
            for (int j = 0; j < bonds.length; ++j) {
                int a1 = mol.indexOf(bonds[j].getAtom1());
                int a2 = mol.indexOf(bonds[j].getAtom2());

                List<Path> path1 = getShortestPaths (i, a1);
                for (Path p1 : path1) {
                    BitSet visited = new BitSet (atoms.length);
                    visited.or(p1.getAtomSet());
                    /*
                     * Unlike the original Horton's algorithm (which was 
                     * devised for minimum basis cycles), we can't
                     * use shortest path here as it fails test2 (i.e.,
                     * it doesn't capture ring size 10).
                     */
                    List<Path> path2 = getPaths (i, a2, visited);
                    for (Path p2 : path2) 
                        if (p1.overlaps(p2, i)) {
                            Ring ring = new Ring (p1, p2, new Path (bonds[j]));
                            if (maxsize <= 2 || ring.size() <= maxsize) {
                                rings.add(ring);
                                if (rings.size() >= maxrings) {
                                    logger.warning("## max number of rings ("
                                                   +maxrings+") reached!");
                                    // stop
                                    return this;
                                }
                            }
                        }
                }
            }

        return this;
    }

    /**
     * Return all rings of sizes <= maxsize
     */
    public Ring[] getRings () { return rings.toArray(new Ring[0]); }
    public int getRingCount () { return rings.size(); }

    // return the number of rings that contains this atom
    public int getRingCount (int atom) { 
        int c = 0;
        for (Ring r : rings)
            if (r.hasAtom(atom)) ++c;
        return c;
    }

    public int getRingCount (MolAtom atom) {
        int idx = mol.indexOf(atom);
        if (idx >= 0) {
            return getRingCount (idx);
        }
        return 0;
    }

    // count the number of ring bonds this atom is connected to
    public int countRingBonds (MolAtom atom) { 
        int c = 0;
        for (int i = 0; i < atom.getBondCount(); ++i) {
            MolBond b = atom.getBond(i);
            int idx = mol.indexOf(b);
            for (Ring r : rings)
                if (r.hasBond(idx)) {
                    ++c;
                    break;
                }
        }
        return c;
    }

    public int getRingCount (MolBond bond) {
        int c = 0;
        int idx = mol.indexOf(bond);
        if (idx >= 0) {
            for (Ring r : rings) 
                if (r.hasBond(idx)) ++c;
        }
        return c;
    }

    /**
     * Return all shortest paths between atoms start and end
     */
    private List<Path> getShortestPaths (int start, int end) {
        List<Path> paths = new ArrayList<Path>();
        dfs1 (paths, new LinkedList<MolBond>(), 
              new BitSet (atoms.length), start, start, end);
        return paths;
    }

    private List<Path> getPaths (int start, int end, BitSet visited) {
        List<Path> paths = new ArrayList<Path>();
        dfs2 (paths, new LinkedList<MolBond>(), visited, start, start, end);
        return paths;
    }
 
    private void dfs1 (List<Path> paths, LinkedList<MolBond> path, 
                       BitSet visited, int start, int a, int end) {
        if (a == end) {
            paths.add(new Path (path));
            return;
        }
        
        visited.set(a);
        MolAtom atom = atoms[a];
        for (int b = 0; b < atom.getBondCount(); ++b) {
            MolBond bond = atom.getBond(b);
            int xa = mol.indexOf(bond.getOtherAtom(atom));

            if ((cost[start][xa] + cost[xa][end] <= cost[start][end])
                && !visited.get(xa)) {
                path.push(bond);
                dfs1 (paths, path, visited, start, xa, end);
                path.pop();
            }
        }
        visited.clear(a);
    }

    private void dfs2 (List<Path> paths, LinkedList<MolBond> path, 
                       BitSet visited, int start, int a, int end) {
        if (a == end) {
            paths.add(new Path (path));
            return;
        }
        
        visited.set(a);
        MolAtom atom = atoms[a];
        for (int b = 0; b < atom.getBondCount(); ++b) {
            MolBond bond = atom.getBond(b);
            int xa = mol.indexOf(bond.getOtherAtom(atom));

            if (!visited.get(xa) && path.size() < maxsize) {
                path.push(bond);
                dfs2 (paths, path, visited, start, xa, end);
                path.pop();
            }
        }
        visited.clear(a);
    }

    static void perception (String mol) throws Exception {
        MolHandler mh = new MolHandler (mol);
        RingPerception rp = new RingPerception ();
        rp.setMaxRings(1024);

        Ring[] rings =  rp.perceive(mh.getMolecule()).getRings();
        System.out.println("..."+rings.length+" rings!");
        for (Ring r : rings) {
            System.out.println(r.size()+": "+r.getAtomSet());
        }
    }

    static void load (InputStream is) throws Exception {
        MolImporter mi = new MolImporter (is);
        long start, time;
        RingPerception rp = new RingPerception ();
        for (Molecule mol = new Molecule (); mi.read(mol); ) {
            System.out.println(">> "+mol.getName());
            start = System.currentTimeMillis();
            Ring[] rings = rp.perceive(mol).getRings();
            time = System.currentTimeMillis() - start;
            System.out.println("..."+rings.length+" ring(s) in "+time+"ms!");
            for (Ring r : rings) {
                System.out.println(String.format("%1$4d: ", r.size())
                                   +" "+r.getAtomSet());
            }
        }
    }

    static void test1 () throws Exception {
        perception (
"\n  Marvin  12201214052D          \n\n"+
" 51 56  0  0  0  0            999 V2000\n"+
"  -35.3479  -10.3420    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n"+
"  -35.3479   -9.5170    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -35.3479   -8.6920    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -34.5229   -9.5170    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -36.1729   -9.5170    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -20.6009  -11.4466    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -20.7472  -10.6347    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -20.1763  -10.0391    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -20.5663   -9.3121    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -20.1316   -8.6109    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -20.5216   -7.8839    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -20.0870   -7.1827    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -19.2624   -7.2084    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -18.8724   -7.9354    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -19.3070   -8.6367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -18.9171   -9.3637    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -19.3517  -10.0649    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -21.3782   -9.4584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -21.4388   -8.6356    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -22.2028   -9.4326    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -21.4900  -10.2758    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -22.2170  -10.6658    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -22.9182  -10.2311    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -23.6453  -10.6211    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -24.3465  -10.1865    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -25.0735  -10.5765    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -25.7747  -10.1419    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -25.7490   -9.3173    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -25.0220   -8.9273    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -25.5323   -8.2791    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -24.4722   -8.3122    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -24.3207   -9.3619    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -26.5017  -10.5319    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -27.2030  -10.0973    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -27.9300  -10.4873    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -28.6312  -10.0527    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -28.6918   -9.2299    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0\n"+
"  -28.0618   -8.6973    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -29.4930   -9.0332    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -29.9276   -9.7345    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -30.7522   -9.7087    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -31.1868  -10.4100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -32.0114  -10.3842    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -32.4014   -9.6572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -31.9668   -8.9560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -31.1422   -8.9817    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -30.7076   -8.2805    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -29.8830   -8.3062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -29.3950  -10.3645    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -29.9597  -10.9660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  -29.1730  -11.1591    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  1  2  1  0  0  0  0\n  2  3  2  0  0  0  0\n  2  4  2  0  0  0  0\n"+
"  2  5  2  0  0  0  0\n  6  7  1  0  0  0  0\n  7  8  1  0  0  0  0\n"+
"  8  9  2  0  0  0  0\n  9 10  1  0  0  0  0\n 10 11  2  0  0  0  0\n"+
" 11 12  1  0  0  0  0\n 12 13  2  0  0  0  0\n 13 14  1  0  0  0  0\n"+
" 14 15  2  0  0  0  0\n 10 15  1  0  0  0  0\n 15 16  1  0  0  0  0\n"+
" 16 17  1  0  0  0  0\n  8 17  1  0  0  0  0\n  9 18  1  0  0  0  0\n"+
" 18 19  1  0  0  0  0\n 18 20  1  0  0  0  0\n 18 21  1  0  0  0  0\n"+
"  7 21  1  0  0  0  0\n 21 22  2  0  0  0  0\n 22 23  1  0  0  0  0\n"+
" 23 24  2  0  0  0  0\n 24 25  1  0  0  0  0\n 25 26  2  0  0  0  0\n"+
" 26 27  1  0  0  0  0\n 27 28  1  0  0  0  0\n 28 29  1  0  0  0  0\n"+
" 29 30  1  0  0  0  0\n 29 31  1  0  0  0  0\n 29 32  1  0  0  0  0\n"+
" 25 32  1  0  0  0  0\n 27 33  2  0  0  0  0\n 33 34  1  0  0  0  0\n"+
" 34 35  2  0  0  0  0\n 35 36  1  0  0  0  0\n 36 37  2  0  0  0  0\n"+
" 37 38  1  0  0  0  0\n 37 39  1  0  0  0  0\n 39 40  2  0  0  0  0\n"+
" 40 41  1  0  0  0  0\n 41 42  2  0  0  0  0\n 42 43  1  0  0  0  0\n"+
" 43 44  2  0  0  0  0\n 44 45  1  0  0  0  0\n 45 46  2  0  0  0  0\n"+
" 41 46  1  0  0  0  0\n 46 47  1  0  0  0  0\n 47 48  2  0  0  0  0\n"+
" 39 48  1  0  0  0  0\n 40 49  1  0  0  0  0\n 36 49  1  0  0  0  0\n"+
" 49 50  1  0  0  0  0\n 49 51  1  0  0  0  0\nM  CHG  2   1  -1  37   1\n"+
"M  END\n"
);
    }

    static void test2 () throws Exception {
        perception (
"\n  Marvin  12201222202D          \n\n"+
" 19 22  0  0  0  0            999 V2000\n"+
"   -8.4384    5.3469    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -7.6453    5.1195    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -7.0767    5.7157    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -6.2812    5.5145    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -6.0654    4.7011    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -5.2508    4.4936    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -4.6092    5.0122    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -4.7375    5.8271    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -5.2192    3.6675    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -4.5429    3.1951    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -6.0023    3.4138    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -6.6165    4.0255    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -6.7475    3.1608    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -7.5168    2.9039    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -7.1875    2.1475    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -8.2413    2.5093    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -7.9740    3.5912    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -8.3353    4.3113    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"   -7.4720    4.2784    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  1  2  1  0  0  0  0\n  2  3  1  0  0  0  0\n  3  4  1  0  0  0  0\n"+
"  4  5  1  0  0  0  0\n  5  6  2  0  0  0  0\n  6  7  1  0  0  0  0\n"+
"  7  8  1  0  0  0  0\n  6  9  1  0  0  0  0\n  9 10  2  0  0  0  0\n"+
"  9 11  1  0  0  0  0\n 11 12  1  0  0  0  0\n  5 12  1  0  0  0  0\n"+
" 12 13  1  0  0  0  0\n 13 14  1  0  0  0  0\n 14 15  1  0  0  0  0\n"+
" 14 16  1  0  0  0  0\n 14 17  1  0  0  0  0\n 17 18  1  0  0  0  0\n"+
" 18 19  1  0  0  0  0\n  2 19  1  0  0  0  0\n 12 19  1  0  0  0  0\n"+
" 17 19  1  0  0  0  0\nM  END\n"
);
    }

    static void test3 () throws Exception {
        perception (
"\n  Marvin  12211211362D          \n\n"+
" 30 38  0  0  0  0            999 V2000\n"+
"    2.0477   -0.6691    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    2.5878   -0.0454    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.3980   -0.2014    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.9381    0.4223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    4.7482    0.2664    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    5.0182   -0.5132    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    5.8284   -0.6691    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    6.3685   -0.0454    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    7.1786   -0.2014    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    7.4486   -0.9809    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    6.9086   -1.6046    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    6.0984   -1.4486    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    5.5583   -2.0723    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    5.8284   -2.8518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    6.6385   -3.0077    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    7.1786   -2.3841    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    5.2883   -3.4755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    4.4782   -3.3196    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.9381   -3.9432    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.1279   -3.7873    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    2.8579   -3.0077    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.3980   -2.3841    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    4.2081   -2.5400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    4.7482   -1.9164    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    4.4782   -1.1368    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.6680   -0.9809    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    3.1279   -1.6046    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    2.3178   -1.4486    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    6.0984    0.7341    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"    5.2883    0.8900    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"+
"  1  2  2  0  0  0  0\n  2  3  1  0  0  0  0\n  3  4  2  0  0  0  0\n"+
"  4  5  1  0  0  0  0\n  5  6  2  0  0  0  0\n  6  7  1  0  0  0  0\n"+
"  7  8  2  0  0  0  0\n  8  9  1  0  0  0  0\n  9 10  2  0  0  0  0\n"+
" 10 11  1  0  0  0  0\n 11 12  2  0  0  0  0\n  7 12  1  0  0  0  0\n"+
" 12 13  1  0  0  0  0\n 13 14  2  0  0  0  0\n 14 15  1  0  0  0  0\n"+
" 15 16  2  0  0  0  0\n 11 16  1  0  0  0  0\n 14 17  1  0  0  0  0\n"+
" 17 18  2  0  0  0  0\n 18 19  1  0  0  0  0\n 19 20  2  0  0  0  0\n"+
" 20 21  1  0  0  0  0\n 21 22  2  0  0  0  0\n 22 23  1  0  0  0  0\n"+
" 18 23  1  0  0  0  0\n 23 24  2  0  0  0  0\n 13 24  1  0  0  0  0\n"+
" 24 25  1  0  0  0  0\n  6 25  1  0  0  0  0\n 25 26  2  0  0  0  0\n"+
"  3 26  1  0  0  0  0\n 26 27  1  0  0  0  0\n 22 27  1  0  0  0  0\n"+
" 27 28  2  0  0  0  0\n  1 28  1  0  0  0  0\n  8 29  1  0  0  0  0\n"+
" 29 30  2  0  0  0  0\n  5 30  1  0  0  0  0\nM  END\n"
);
    }

    public static void main (String[] argv) throws Exception {
        //test1 ();
        //test2 ();
        //test3 ();
        if (argv.length == 0) {
            logger.info("## reading from STDIN...");
            load (System.in);
        }
        else {
            for (String a : argv) {
                load (new FileInputStream (new File (a)));
            }
        }
    }
}