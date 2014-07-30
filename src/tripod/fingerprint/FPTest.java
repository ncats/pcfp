package tripod.fingerprint;

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

import java.io.*;
import java.util.*;
import java.util.zip.*;

import chemaxon.struc.Molecule;
import chemaxon.formats.MolImporter;
import chemaxon.util.MolHandler;

public class FPTest {
    int maxdiff = 0, total = 0;
    long totaldif = 0;
    Molecule maxmol;
    PCFP pcfp = new PCFP ();
    int[] fn = new int[PCFP.FP_SIZE];
    int[] fp = new int[PCFP.FP_SIZE];

    FPTest () {
    }
    
    void evaluate (PrintStream ps, Molecule mol) throws Exception {
        mol.hydrogenize(false);
        long start = System.currentTimeMillis();
        pcfp.setMolecule(mol);
        long end = System.currentTimeMillis();

        String enc = mol.getProperty("PUBCHEM_CACTVS_SUBSKEYS");
        if (enc == null) {
            return;
        }
        ++total;

        PCFP fp = PCFP.decode(enc);
        int diff = pcfp.hamming(fp);
        if (diff == 0) {
            String b64 = pcfp.encode();
            if (!b64.equals(enc)) {
                throw new IllegalStateException
                    ("Expecting "+enc+" but got "+b64);
            }
            return;
        }

        // diff != 0
        BitSet bits2 = fp.toBits();
        bits2.andNot(pcfp.toBits());

        // for now only check if we don't have the bits that we should
        // have (i.e., false negatives); not too concern about false
        // positives.
        if (bits2.cardinality() > 0) {
            BitSet bits = pcfp.toBits();
            ps.println(mol.getName()+" "+bits.cardinality()+" "+bits);
            bits.andNot(fp.toBits());
            ps.print("## time = "+(end-start)+"ms");
            ps.println("  +"+bits+" -"+bits2);
            for (int i = bits2.nextSetBit(0); 
                 i >= 0; i = bits2.nextSetBit(i+1)) 
                ++this.fn[i];
            for (int i = bits.nextSetBit(0); i>=0; i = bits.nextSetBit(i+1))
                ++this.fp[i];
            
            ps.println("## Hamming distance = " + diff);
            ps.println("-- PubChem's original bits --");
            bits2 = fp.toBits();
            ps.println(mol.getName()+" "+bits2.cardinality()+" "+bits2);
            
            ps.println(">> "+mol.toFormat("smiles:q") +"\t"+mol.getName());
            ps.println();
        }

        if (diff > maxdiff || maxmol == null) {
            maxdiff = diff;
            maxmol = mol;
        }
        totaldif += diff;
    }

    void summary () {
        if (maxmol != null) {
            double err = (double)totaldif/total;
            System.out.print("** false-negative:");
            for (int i = 0; i < fn.length; ++i) 
                if (fn[i] > 0) System.out.print(" "+i+":"+fn[i]);
            System.out.println();
            System.out.print("** false-positive:");
            for (int i = 0; i < fp.length; ++i) 
                if (fp[i] > 0) System.out.print(" "+i+":"+fp[i]);
            System.out.println();
            System.out.println("** total diff: "+totaldif);
            System.out.println("** total molecules: "+total);
            System.out.println("** average err per molecule: " + err);
            System.out.println("** maxdiff: "+maxdiff);
            System.out.println
                ("** worst offender: "
                 +maxmol.toFormat("smiles:q")+"\t"+maxmol.getName());
        }
    }

    void generate (PrintStream ps, String id, Molecule mol) {
        pcfp.setMolecule(mol);
        ps.print(mol.getName());
        BitSet bits = pcfp.toBits();
        for (int i = bits.nextSetBit(0); i >= 0; i = bits.nextSetBit(i+1)) {
            ps.print(" "+i);
        }
        ps.println();
    }

    void shutdown () {
        pcfp.shutdown();
    }

    public static void main (String argv[]) throws Exception {
	if (argv.length == 0) {
	    System.out.println("FPTest FILES...");
	    System.exit(1);
	}

        FPTest fptest = new FPTest ();

        int id = 0;
	for (int i = 0; i < argv.length; ++i) {
	    MolImporter importer = null;
	    try {
		importer = new MolImporter 
		    (new GZIPInputStream (new FileInputStream (argv[i])));
	    }
	    catch (Exception ex) {
                try {
                    importer = new MolImporter (argv[i]);
                }
                catch (Exception exx) {
                    fptest.evaluate
                        (System.out, new MolHandler (argv[i]).getMolecule());
                }
	    }

            if (importer != null) {
                for (Molecule mol = new Molecule (); importer.read(mol);) {
                    //fptest.evaluate(System.out, mol);
                    fptest.generate(System.out, String.valueOf(++id), mol);
                }
                importer.close();
            }
	}
        fptest.summary();
        fptest.shutdown();
    }
}
