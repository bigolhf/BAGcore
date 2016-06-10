/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core;

import java.util.ArrayList;

/**
 * a simple sequence object that contains some basic functionality
 * and stores basic sequence features
 * 
 * @author weibo
 */
public class SimpleSeq {
    
    
    private String          id="";
    private String          seq="";
    private int             length=0;
    private int             start=0;
    private int             end=0;
    private String          chromosome="";
    private String          strand;
    private String          name="";
    private String          accessionNumber;


    /**
     * Empty Class Constructor
     */
    public SimpleSeq(){

    }
    
    public SimpleSeq(SimpleSeq nSimpleSeq){
        this.chromosome         = nSimpleSeq.chromosome;
        this.accessionNumber    = nSimpleSeq.accessionNumber;
        this.strand             = nSimpleSeq.strand;
        this.start              = nSimpleSeq.start;
        this.end                = nSimpleSeq.end;
        this.length             = nSimpleSeq.length;
        this.name               = nSimpleSeq.name;
        this.id                 = nSimpleSeq.id;
        this.seq                = nSimpleSeq.seq;
    }
    
    
    
    /**
     * Constructor 
     * 
     * @param id
     * @param seq 
     */
    public SimpleSeq(String id,String seq){
        this.name   = id;
        this.id     = id;
        this.seq    = seq;
        this.length = seq.length();
        
    }
    
    
    
    
    /**
     * summarize the sequence properties
     * 
     * @return 
     */
    @Override
    public String toString(){
        String str = "ID:" + this.getId() + "\t"
                    + "Name::" + this.getName() + "\t"
                    + "Start:" + this.getStart() + "\t"
                    + "End:" + this.getEnd() + "\t"
                    + "Len:" + this.getLength() + "\t"
          + "\n";
        return str;
    }


    /**
     * return subsequence at specified position
     * 
     * @param start
     * @param stop
     * @return 
     */
    public String getSequence(int start, int stop){
        if(start>=0 && stop < seq.length())
            return seq.substring(start, stop+1);

        return null;
    }

    
    /**
     * convert RNA sequence to DNA
     * 
     * @param seq
     * @return 
     */
    public static String rna2dna(String rnaseq){
        return rnaseq.replace("u", "t").replace("U", "T");
    }
    
    
    /**
     * convert DNA sequence to RNA
     * 
     * @param dnaseq
     * @return 
     */
    public static String dna2rna(String dnaseq){
        return dnaseq.replace("t", "u").replace("T", "U");
    }
    
    
    /**
     * Reverse complement RNA/DNA sequence
     * if sequence contains 'U' or 'u' it is assumed the sequence is RNA
     * 
     * @param seqIn
     * @return 
     */
    public static String complement(String seqIn){
        
        StringBuilder Complement = new StringBuilder();
        char [] strReversed = new StringBuilder(seqIn).reverse().toString().toCharArray();

        for (char nt: strReversed) {
            switch (nt){
                case 'a':
                    if(seqIn.toLowerCase().contains("u"))
                        Complement.append("t");
                    else
                        Complement.append("u");
                    break;
                    
                case 'A':
                    if(seqIn.toLowerCase().contains("u"))
                        Complement.append("T");
                    else
                        Complement.append("U");
                    break;
                    
                case 'c':
                    Complement.append("g");
                    break;
                    
                case 'C':
                    Complement.append("G");
                    break;
                    
                case 'g':
                    Complement.append("c");
                    break;
                    
                case 'G':
                    Complement.append("C");
                    break;
                    
                case 't':
                    Complement.append("a");
                    break;
                    
                case 'T':
                    Complement.append("A");
                    break;
                    
                case 'u':
                    Complement.append("a");
                    break;
                    
                case 'U':
                    Complement.append("A");
                    break;
                    
                default:
                    Complement.append("N");
                    break;
                    
            }

        }
        return Complement.toString();
    }


    /**
     * count occurrence of a specific nucleotide within the specified sequence
     * 
     * @param qSeq
     * @param nt
     * @return 
     */
    public static int NTcount(String qSeq, char nt){
        char[] cs=qSeq.toCharArray();
        int n=0;
        for(char c:cs){
            if(c==nt){
                n++;
            }
        }
        return n;
    }
    
    
    
    
    
    /**
     * split the sequence into fragments according to specified step and window size
     * 
     * @param step
     * @param window
     * @return 
     */
    public ArrayList<SimpleSeq> splitSequence(int step, int window){
        
        ArrayList<SimpleSeq> fragmentList = new ArrayList<>();
        length=seq.length();
        int n=(length-window)/step+1;
        if(n<1) n=1;

        int fragStart=0;
        int fragEnd; 
        for(int i=0;i<n;i++){
            
            if(fragStart>=length) break;
            fragEnd = fragStart + window;
            if(fragEnd>length) fragEnd=length;
            String fragID = name + "_" + (fragStart+1) + "-"+fragEnd;
            String subseq = seq.substring(fragStart,fragEnd);

            SimpleSeq frag = new SimpleSeq(fragID,subseq);
            frag.setAbsStartInQuerySeq(fragStart+1);// count from 1
            frag.setAbsEndInQuerySeq(fragEnd); //count from 1
            frag.setName(fragID);
            fragmentList.add(frag);
            fragStart+= step;
            
        }
        return fragmentList;
    }
    
    
    
    /**
     * count occurrence of a specific nucleotide within the defined sequence
     * 
     * @param nt
     * @return 
     */
    public int NTcount(char nt){
        return NTcount(seq, nt);
    }

    /**
     * @return the id
     */
    public String getId() {
        return id;
    }

    /**
     * @param id 
     */
    public void setID(String id) {
        this.id = id;
    }

    /**
     * @return the seq
     */
    public String getSeq() {
        return seq;
    }

    /**
     * @param seq 
     */
    public void setSeq(String seq) {
        this.seq = seq;
    }

    /**
     * @return length
     */
    public int getLength() {
        return length;
    }

    /**
     * @param length 
     */
    public void setLength(int length) {
        this.length = length;
    }

    /**
     * @return start
     */
    public int getStart() {
        return start;
    }

    /**
     * @param start
     */
    public void setAbsStartInQuerySeq(int start) {
        this.start = start;
    }

    /**
     * @return end
     */
    public int getEnd() {
        return end;
    }

    /**
     * @param end
     */
    public void setAbsEndInQuerySeq(int end) {
        this.end = end;
    }

    /**
     * @return name
     */
    public String getName() {
        return name;
    }

    /**
     * @param name 
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @return the accessionNumber
     */
    public String getAccessionNumber() {
        return accessionNumber;
    }

    /**
     * @param accessionNumber the accessionNumber to set
     */
    public void setAccessionNumber(String accessionNumber) {
        this.accessionNumber = accessionNumber;
    }

    /**
     * @return the chromosome
     */
    public String getChromosome() {
        return chromosome;
    }

    /**
     * @param chromosome the chromosome to set
     */
    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    /**
     * @return the strand
     */
    public String getStrand() {
        return strand;
    }

    /**
     * @param strand the strand to set
     */
    public void setStrand(String strand) {
        this.strand = strand;
    }

    


}