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


public class SearchParams implements java.io.Serializable {
    private static final long serialVersionUID = 0x07424ca4cf6c8e9al;

    public enum Type {
        Substructure,
            Superstructure,
            Similarity,
            Exact
            };


    // hard limit on screen results
    protected static final int DEFAULT_MATCH_LIMIT = 10000;
    // minimum fingerprint density required to proceed
    protected static final double DEFAULT_MIN_DENSITY = 0.005;
    // milliseconds
    protected static final long DEFAULT_SEARCH_TIMEOUT = 1000*60;
    protected static final double DEFAULT_SIMILARITY = 0.7;

    protected int matchLimit = DEFAULT_MATCH_LIMIT;
    protected double minDensity = DEFAULT_MIN_DENSITY;
    protected long searchTimeout = DEFAULT_SEARCH_TIMEOUT;
    protected double similarity = DEFAULT_SIMILARITY;
    protected boolean aromatize = true;
    protected Type type;

    public SearchParams () {
        this (Type.Substructure);
    }

    public SearchParams (Type type) {
        this.type = type;
    }

    public static SearchParams substructure () { 
        return new SearchParams (Type.Substructure); 
    }
    public static SearchParams superstructure () {
        return new SearchParams (Type.Superstructure);
    }
    public static SearchParams similarity () {
        return new SearchParams (Type.Similarity);
    }
    public static SearchParams similarity (double similarity) {
        return new SearchParams (Type.Similarity).setSimilarity(similarity);
    }
    public static SearchParams exact () {
        return new SearchParams (Type.Exact);
    }

    public SearchParams setType (Type type) { 
        this.type = type;
        return this;
    }
    public Type getType () { return type; }

    public void setAromatize (boolean aromatize) { 
        this.aromatize = aromatize;
    }
    public boolean getAromatize () { return aromatize; }

    public SearchParams setMatchLimit (int limit) {
	this.matchLimit = limit;
        return this;
    }
    public int getMatchLimit () { return matchLimit; }

    public SearchParams setMinDensity (double density) {
	this.minDensity = density;
        return this;
    }
    public double getMinDensity () { return minDensity; }

    public SearchParams setTimeout (long timeout) { // in miliseconds
        searchTimeout = timeout;
        return this;
    }
    public long getTimeout () { return searchTimeout; }

    public SearchParams setSimilarity (double similarity) {
        this.similarity = similarity;
        return this;
    }
    public double getSimilarity () { return similarity; }

    public String toString () {
        return "{type="+type+",limit="+matchLimit+",density="+minDensity
            +",timeout="+searchTimeout+",similarity="+similarity+"}";
    }
}
