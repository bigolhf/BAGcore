/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core;

import java.util.ArrayList;
import java.util.HashMap;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 * stores an miRNA entry, while this may seem a cumbersome implementation, the most
 * miRNAs found in a mammalian is less than 3000, so we shouldn't be too worried 
 * about overheads
 * 
 * This needs to be made comparable so the entries stored in a list can be sorted
 * to speed things up
 * 
 * @author sr
 */
public class MiRNAFeature extends SimpleRNASequence{

    static Logger               logger                      = LogManager.getLogger();
    
    private static final int        QUALIFIER_START_COL     =   21;
    private static final String     QUALIFIER_ACCESSION     =   "/accession";
    private static final String     QUALIFIER_PRODUCT       =   "/product";
    private static final String     QUALIFIER_EVIDENCE      =   "/evidence";
    private static final String     QUALIFIER_EXPERIMENT    =   "/experiment";
    private String              mimatID;
    private String              note;
    private String              name;
    private String              parent;
    private String              chromosome;
    private int                 startPos;
    private int                 endPos;
    private String              strand;
    private String              seq;
    private String              isomiRString;
    private HashMap             featureSet = new HashMap();     // stores additional characteristics of this miRNA
    
    private String              evidence="";
    private String              references="";                     // string of ";" delimited pubmed IDs
    
    
    private String              host;
    
    private int                 largestInternalBulge;
    private int                 numOfBulges;
    private int                 numOfUnpairedBases;
    private double              fractOfUnpairedBases;

    private char                firstBase;
    private double              stability;
    private char                dangleBaseOne;
    private char                dangleBaseTwo;

    private int                 miStart;
    private int                 miEnd;

    
    // for storing isomiR information
    //(name + ";" + start + ";" + cigar + ";" + md + ";" + seq + "\t")
    
    private static final int    NAMECOL    = 0;
    private static final int    STARTCOL   = 1;
    private static final int    CIGARCOL   = 2;
    private static final int    MDCOL      = 3;
    private static final int    SEQCOL     = 4;
    
    

    public MiRNAFeature(){
        
    }
    
    /**
     * Constructor specifying name, location and sequence of an miRNA
     * (ideal for predicted miRNAs)
     * 
     * @param cName     String      name
     * @param cChrom    String      chromosome
     * @param cStart    int         start position
     * @param cEnd      int         end position
     * @param cStrand   String      Strand
     * 
     */
    public MiRNAFeature(String cName, String cChrom, int cStart, int cEnd, String cStrand){
        this(cName, cChrom, cStart, cEnd, cStrand, "", "", "");
    }
    
    
    
    /**
     * Constructor specifying name, location and sequence of an miRNA
     * (ideal for miRBase entry miRNAs)
     * 
     * @param cName     String      name
     * @param cChrom    String      chromosome
     * @param cStart    int         start position
     * @param cEnd      int         end position
     * @param cStrand   String      Strand
     * @param cMIMAT    String      MIMAT (miRBase) ID
     * @param cParentID String      MI (miRBase) Parent ID
     * @param cSeq      String      sequence
     * 
     */
    public MiRNAFeature(String cName, String cChrom, int cStart, int cEnd, String cStrand, String cMIMAT, String cParentID, String cSeq){
        name = cName;
        chromosome = cChrom;
        startPos = cStart;
        endPos = cEnd;
        strand = cStrand;
        parent = cParentID;
        mimatID =cMIMAT;
        seq = cSeq;
        isomiRString = "";
        
    }
    
    
    public MiRNAFeature(MiRNAFeature m){
        
        name            = m.name;
        chromosome      = m.chromosome;
        startPos        = m.startPos;
        endPos          = m.endPos;
        strand          = m.strand;
        parent          = m.parent;
        mimatID         = m.mimatID;
        seq             = m.seq;
        isomiRString    = m.isomiRString;
        
    }
    
    
    
    
    /**
     * Add feature to set
     * 
     * @param key
     * @param value
     * @return 
     */
    public int addFeature(String key, String value){
        featureSet.put(key, value);
        return featureSet.size();
    }
    
    
    
    
    
    /**
     * checks whether Chromosome strings are the same, while attempting
     * to allow for the presence or absence of a variation on the 'Chr' 
     * prefix
     * 
     * @param queryChr
     * @return 
     */
    public Boolean chromosomeMatch(String queryChr){
        
        return MiRNAFeature.removeChromosomePrefix(chromosome).equals(MiRNAFeature.removeChromosomePrefix(queryChr));
        
    }
    
    
    
    
    /**
     * attempt to remove any prefix of the form 'Chr' from the chromosome string
     * 
     * @param chrString
     * @return 
     */
    public static String removeChromosomePrefix(String chrString){
        
        if(chrString.contains("chr")){
            chrString = chrString.replace("chr", "");
        }
        else{
            if(chrString.contains("CHR")){
                chrString = chrString.replace("CHR", "");
            }
            else{
                if(chrString.contains("Chr")){
                    chrString = chrString.replace("Chr", "");
                }
            }
            
        }
        return chrString;
        
    }
        
    
    /**
     * characterize this miRNA
     * some of the features can only be characterized from the parent pri-miRNA
     * 
     * @return 
     */
    public HashMap characterize(){
        HashMap characteristics = new HashMap();
        
        return characteristics;
    }
    
    
    
    
    
    
    
    
    
    
    
    /**
     * add information to define an isomiR for this entry
     * 
     * @param name
     * @param start
     * @param cigar
     * @param md 
     * @param seq
     * 
     */
    public void addIsomiR(String name, int start, String cigar, String md, String seq){

        isomiRString = isomiRString.concat(name + ";" + start + ";" + cigar + ";" + md + ";" + seq + "\t");
        
    }
    
    /**
     * delete the isomiR string
     * 
     */
    public void removeIsomiRs(){
        isomiRString = "";
    }

    
    /**
     * write out isomiRs.
     * 
     * only report reads that are above a baseline, defined in terms of the fraction 
     * of the total number of reads for the miRNA. e.g. if there are 100 reads, and 
     * baseline is 5, then only isomiRs with more than 5 reads will be reported
     * 
     * @param baselinePercent : int     only report isomiRs with reads that
     * @param minCounts       : int     total counts for isomiR must be greater
     *                                  than this value
     * @return 
     */
    public String reportIsomiRs(int baselinePercent, int minCounts){
        String reportStr = this.getName() + ":\t[" + this.getChromosome() + "]\t" 
                + this.getStartPos() + "\t" + this.getEndPos() + this.getSequence() + "\n";
        String [] isomiRs = isomiRString.split("\t");
        
        int totalCounts = this.getTotalCounts();
        if (totalCounts < minCounts) return "";
        reportStr = reportStr.concat("Total Counts = " + totalCounts + "\n");
        
        String isomiRStr = "";
        for(String isomiR: isomiRs){
            String[] values = isomiR.split(";");
            if(Double.parseDouble(isomiR.split(";")[NAMECOL].split("-")[1]) / (double) totalCounts > (double)baselinePercent/100.0){
                isomiRStr = isomiRStr.concat("name: " + values[0] + "\t"
                            + "start: " + values[STARTCOL] + "\t"
                            + "cigar: " + values[CIGARCOL] + "\t"
                            + "MD: " + values[MDCOL] + "\t" 
                            + "SQ: " + values[SEQCOL] + "\n"
                );
            }
        }
        
        if (isomiRStr.equals("")) return "";
        return reportStr.concat(isomiRStr);
    }
    
    
    
    
    /**
     * report the isomiR in a manner that is visually appealing
     * 
     * @param baselinePercent
     * @param minCounts
     * @return 
     */
    public String prettyReportIsomiRs(int baselinePercent, int minCounts){
        
        String reportStr = this.getName() + "|" + this.getMimatID() + " :\tchr" + this.getChromosome() + "\t" 
                + this.getStartPos() + " --> " + this.getEndPos() + " (" + this.getStrand() + ") : " + this.getSequence() + "\n";
        
        int totalCounts = this.getTotalCounts();
        reportStr = reportStr.concat("Total Counts = " + totalCounts + "\n");
        logger.debug(reportStr);
        if (totalCounts < minCounts) return "";
        
    
        String [] isomiRs = isomiRString.split("\t");
        
        int longestName = this.getName().length();
        int longestSeq = this.getSequence().length();
        int longestCounts = 0;
        int longestMD = 0;
        int minStart = this.getIsomiRMinStart();
        int maxStop = this.getIsomiRMaxStop();
        
        for(String isomiR: isomiRs){
            
            if(Double.parseDouble(isomiR.split(";")[NAMECOL].split("-")[1]) / (double) totalCounts > (double)baselinePercent/100.0){
                if(isomiR.split(";")[NAMECOL].length() > longestName)
                    longestName = isomiR.split(";")[NAMECOL].length();

                if(isomiR.split(";")[NAMECOL].split("-")[1].length() > longestCounts)
                    longestCounts = isomiR.split(";")[NAMECOL].split("-")[1].length();

                if(isomiR.split(";")[SEQCOL].length() > longestSeq)
                    longestSeq = isomiR.split(";")[SEQCOL].length();
                if(this.getStrand().equals("+")){
                    if(isomiR.split(";")[SEQCOL].equals(this.getSequence().replace("U", "T")))
                        longestName++;
                }
                else{
                    if(isomiR.split(";")[SEQCOL].equals(SimpleSeq.complement(this.getSequence()).replace("U", "T")))
                        longestName++;                    
                }
                
                    

                if(isomiR.split(";")[MDCOL].length() > longestMD)
                    longestMD = isomiR.split(";")[MDCOL].length();
                
                
            }          

        }    
        
        int leftMargin = 10;
        int ColMargin = 5;
        
        String isoString = "";
        for(String isomiR: isomiRs){
            
            if(Double.parseDouble(isomiR.split(";")[NAMECOL].split("-")[1]) / (double) totalCounts > (double)baselinePercent/100.0){

                isoString           = isoString.concat(StringUtils.repeat(" ", leftMargin));
                String isoName      = isomiR.split(";")[NAMECOL];
                if(this.getStrand().equals("+")){
                    if(isomiR.split(";")[SEQCOL].equals(this.getSequence().replace("U", "T")))
                        isoName = "*".concat(isoName);
                }
                else{
                    if(isomiR.split(";")[SEQCOL].equals(SimpleSeq.complement(this.getSequence()).replace("U", "T")))
                        isoName = "*".concat(isoName);
                }
                isoString           = isoString.concat(StringUtils.repeat(" ", longestName - isoName.length()) + isoName);
                
                int isoStart        = Integer.parseInt(isomiR.split(";")[STARTCOL]);
                String isoSeq       = isomiR.split(";")[SEQCOL];
                isoString           = isoString.concat(StringUtils.repeat(" ", ColMargin + (isoStart - this.getIsomiRMinStart())));
                isoString           = isoString.concat(isoSeq);
                isoString           = isoString.concat(StringUtils.repeat(" ", (this.getIsomiRMaxStop()-(isoStart + isoSeq.length())) + ColMargin));
                
                String isoMD        = isomiR.split(";")[MDCOL];
                isoString           = isoString.concat(StringUtils.repeat(" ", ColMargin) + isoMD + StringUtils.repeat(" ", (longestMD - isoMD.length()) + ColMargin));
                
                String isoCounts    = isomiR.split(";")[NAMECOL].split("-")[1];
                isoString           = isoString.concat(StringUtils.repeat(" ", ColMargin + (longestCounts - isoCounts.length())) + isoCounts + StringUtils.repeat(" ", ColMargin) + "\n");
            
            }
            
        }
        
        if (isoString.equals("")) return "";
        
        return reportStr.concat(isoString + "\n\n\n");

    }
    
    
    
    
        /**
     * an isomiR can be classified as 5´ modification , 3´ modification or polymorphic
     * 
     * @param baselinePercent
     * @param minCounts 
     * 
     * @return ArrayList    : list of isomiR points
     */
    public ArrayList characterizeIsomiRs(int baselinePercent){
        
        ArrayList isomiRPts = new ArrayList<>();
        // we can identify 5´ modification from start position
        int totalCounts = this.getTotalCounts();
        String [] isomiRStrings = isomiRString.split("\t");
        String isoString = "";
        logger.debug(this.name);

        for(String isomiR: isomiRStrings){

            
            if(Double.parseDouble(isomiR.split(";")[NAMECOL].split("-")[1]) / (double) totalCounts > (double)baselinePercent/100.0){
                logger.debug(isomiR.split(";")[NAMECOL] +  this.getStrand() +  "\n" + " : " + Double.parseDouble(isomiR.split(";")[NAMECOL].split("-")[1]) / (double) totalCounts);


                //logger.info(isomiR);
                int isoStart        = Integer.parseInt(isomiR.split(";")[STARTCOL]);
                String isoSeq       = isomiR.split(";")[SEQCOL];
                int isoEnd          = isoStart + isoSeq.length() - 1;

                int noOf5pSteps = 0;
                int noOf3pSteps = 0;
                int noOfPolySteps = 0;                

                int wStart = 0;
                int wEnd = 0;
                int iStart = 0;
                int iEnd = 0;

                
                if(this.getStrand().equals("+")){
                    logger.debug(this.seq + "\n" + isomiR.split(";")[SEQCOL]);
                    if(isoStart != this.startPos){                   
                        noOf5pSteps = isoStart - this.startPos;                   
                    }

                    if(isoEnd != this.endPos){
                        noOf3pSteps = isoEnd - this.getEndPos();
                    }
                    logger.debug(noOf5pSteps + ", " + noOf3pSteps);
                    logger.debug(isoStart + ", " + isoEnd + ", " + (this.startPos- isoStart) + ", " + seq.length());

                    if(this.startPos >= isoStart ){
                        wStart = 0;
                        iStart = this.startPos- isoStart;
                    }
                    else{
                        wStart = isoStart - this.startPos;
                        iStart = 0;
                    }

                    if(this.endPos < isoEnd){
                        wEnd = seq.length(); 
                        iEnd = isoSeq.length() - (isoEnd - this.endPos);
                    }
                    else{
                        wEnd = seq.length() - (this.endPos - isoEnd);
                        iEnd = isoSeq.length();
                    }

                    if(seq.substring(wStart, wEnd).replace("U", "T").equals(isoSeq.substring(iStart, iEnd))==false){
                        for(int b=0; b<wEnd-wStart; b++){
                            if(seq.charAt(wStart + b) != isoSeq.charAt(iStart + b)) noOfPolySteps++;
                        }
                    }
                }
                else{
                    logger.debug(this.seq + "\t" + startPos + "\t" + endPos + "\n");
                    logger.debug(SimpleSeq.complement(isomiR.split(";")[SEQCOL]) + "\t" + isoStart + "\t" + isoEnd + "\n");
                    
                    if(isoStart != this.startPos){                   
                        noOf5pSteps = isoStart - this.startPos;                   
                    }

                    if(isoEnd != this.endPos){
                        noOf3pSteps = isoEnd - this.getEndPos();
                    }
                    logger.debug(noOf5pSteps + ", " + noOf3pSteps);
                    logger.debug(isoStart + ", " + isoEnd + ", " + (this.startPos- isoStart) + ", " + seq.length());

                    if(this.startPos >= isoStart ){
                        wStart = 0;
                        iStart = this.startPos- isoStart;
                    }
                    else{
                        wStart = isoStart - this.startPos;
                        iStart = 0;
                    }

                    if(this.endPos < isoEnd){
                        wEnd = seq.length(); 
                        iEnd = isoSeq.length() - (isoEnd - this.endPos);
                    }
                    else{
                        wEnd = seq.length() - (this.endPos - isoEnd);
                        iEnd = isoSeq.length();
                    }

                    if(seq.substring(wStart, wEnd).replace("U", "T").equals(SimpleSeq.complement(isoSeq.substring(iStart, iEnd)))==false){
                        String complementWTSeq = SimpleSeq.complement(seq);
                        for(int b=0; b<wEnd-wStart; b++){
                            if(complementWTSeq.charAt(wStart + b) != isoSeq.charAt(iStart + b)) noOfPolySteps++;
                        }
                    }                 
                    
                }


                String isoCounts    = isomiR.split(";")[NAMECOL].split("-")[1];
                
                HashMap isomiRPt = new HashMap();
                isomiRPt.put("5p", noOf5pSteps);
                isomiRPt.put("3p", noOf3pSteps);
                isomiRPt.put("poly", noOfPolySteps);
                isomiRPt.put("fraction", Double.parseDouble(isoCounts)/(double)totalCounts);
                
                isomiRPts.add(isomiRPt);

            }

        }
        
        return isomiRPts;

    }
    

    
    
    
    /**
     * sum counts for this isomiR
     * 
     * @return 
     */
    public int getTotalCounts(){
        
        int totalCounts = 0;
        String [] isomiRs = isomiRString.split("\t");

        for(String isomiR: isomiRs){            
            totalCounts += Integer.parseInt(isomiR.split(";")[NAMECOL].split("-")[1]);
        }
        return totalCounts;
    }
    
    
    /**
     * find the smallest start position within the isomiRs
     * 
     * @return 
     */
    public int getIsomiRMinStart(){
        
        int isomiRMinStart = 1000000000;
        String [] isomiRs = isomiRString.split("\t");

        for(String isomiR: isomiRs){         
            if(Integer.parseInt(isomiR.split(";")[STARTCOL]) < isomiRMinStart) 
                isomiRMinStart = Integer.parseInt(isomiR.split(";")[STARTCOL]);
        }
        return isomiRMinStart;
        
    }
    
    
    
    
    /**
     * find the smallest start position within the isomiRs
     * 
     * @return 
     */
    public int getIsomiRMaxStop(){
        
        int isomiRMaxStop = -1;
        String [] isomiRs = isomiRString.split("\t");

        for(String isomiR: isomiRs){         
            if(Integer.parseInt(isomiR.split(";")[STARTCOL]) + isomiR.split(";")[SEQCOL].length() > isomiRMaxStop) 
                isomiRMaxStop = Integer.parseInt(isomiR.split(";")[STARTCOL]) + isomiR.split(";")[SEQCOL].length();
        }
        return isomiRMaxStop;
        
    }
    
    
    
    /**
     * this will only parse miRNA entries defined in EMBL format in the miRNA.dat
     * file released by miRBase. i.e., only a subset of qualifies are recognized
     * e.g. 
     * 
     *  FT   miRNA           17..38
     *  FT                   /accession="MIMAT0000001"
     *  FT                   /product="cel-let-7-5p"
     *  FT                   /evidence=experimental
     *  FT                   /experiment="cloned [1-3], Northern [1], PCR [4], 454 [5],
     *  FT                   Illumina [6], CLIPseq [7]"
     * 
     * @param emblLines
     */
    public void parseMiRBaseMiRNAlines(ArrayList<String> emblLines){
        String lastQualifier = "";
        for(String emblLine:emblLines){
            // header line
            if(emblLine.split("\\s+")[1].trim().equals("miRNA")){
                this.miStart= Integer.parseInt(emblLine.substring(QUALIFIER_START_COL, emblLine.indexOf("..")));
                this.miEnd=Integer.parseInt(emblLine.substring(emblLine.indexOf("..")+2).trim());
                continue;
            }
            
            // key / value entries
            String qualifier;
            if(emblLine.contains("="))
                qualifier = emblLine.substring(QUALIFIER_START_COL, emblLine.indexOf("="));
            else
                qualifier="";
            
            switch(qualifier){
                case QUALIFIER_ACCESSION:
                    this.mimatID = emblLine.substring(emblLine.indexOf("=")+1).replaceAll("\"", "");
                    lastQualifier = QUALIFIER_ACCESSION;
                    break;
                    
                case QUALIFIER_PRODUCT:
                    this.name = emblLine.substring(emblLine.indexOf("=")+1).replaceAll("\"", "");
                    lastQualifier = QUALIFIER_PRODUCT;
                    break;
                    
                case QUALIFIER_EVIDENCE:
                    this.evidence = emblLine.substring(emblLine.indexOf("=")+1).replaceAll("\"", "");
                    lastQualifier = QUALIFIER_EVIDENCE;
                    break;
                    
                case QUALIFIER_EXPERIMENT:
                    this.references = emblLine.substring(emblLine.indexOf("=")+1).replaceAll("\"", "");
                    lastQualifier = QUALIFIER_EXPERIMENT;
                    break;
                    
                default:
                    concatKeyValue(emblLine, lastQualifier);
            }
        }
    }
    
    
    /**
     * this handles continuation lines in the key/value entry. It should only be 
     * called from @see parseMiRBaseMiRNAlines(ArrayList<String> emblLines). it 
     * is only included in a separate method for readability
     * 
     * @param emblLine
     * @param lastQualifier
     * @return 
     */
    private String concatKeyValue(String emblLine, String lastQualifier){
        
        switch(lastQualifier){
            case QUALIFIER_ACCESSION:
                this.mimatID = mimatID.concat(emblLine.substring(QUALIFIER_START_COL).replaceAll("\"", "") + " ");
                lastQualifier = QUALIFIER_ACCESSION;
                break;

            case QUALIFIER_PRODUCT:
                this.name = name.concat(emblLine.substring(QUALIFIER_START_COL).replaceAll("\"", "") + " ");
                lastQualifier = QUALIFIER_ACCESSION;
                break;

            case QUALIFIER_EVIDENCE:
                this.evidence = evidence.concat(emblLine.substring(QUALIFIER_START_COL).replaceAll("\"", "") + " ");
                lastQualifier = QUALIFIER_ACCESSION;
                break;

            case QUALIFIER_EXPERIMENT:
                this.references = references.concat(emblLine.substring(QUALIFIER_START_COL).replaceAll("\"", "") + " ");
                lastQualifier = QUALIFIER_ACCESSION;
                break;

            default:
      
        }                    
      
        return lastQualifier;
    }
    
    
    
    
    /**
     * 
     * @param miRFeat
     * @return
     */
    @Deprecated
    public int compareTo(MiRNAFeature miRFeat) {

        int thisChr = -1;
        int queryChr = -1;
        if(this.getChromosome().contains("chr")){
            thisChr = Integer.parseInt(this.getChromosome().replace("chr", ""));
            queryChr = Integer.parseInt(miRFeat.getChromosome().replace("chr", ""));
        }
        
        if(this.getChromosome().contains("CHR")){
            thisChr = Integer.parseInt(this.getChromosome().replace("CHR", ""));
            queryChr = Integer.parseInt(miRFeat.getChromosome().replace("CHR", ""));
        }
        
        if(this.getChromosome().contains("Chr")){
            thisChr = Integer.parseInt(this.getChromosome().replace("Chr", ""));
            queryChr = Integer.parseInt(miRFeat.getChromosome().replace("Chr", ""));
        }
        
     return thisChr - queryChr;

     }
   
    
    
    
    /**
     * base equality on the mimatID, which should be unique
     * 
     * @param qObject
     * @return 
     */
    @Override
    public boolean equals(Object qObject){
        if (qObject != null && qObject instanceof MiRNAFeature)
        {
            return (this.mimatID.equals(((MiRNAFeature) qObject).mimatID));
        }
        return false;

    }
    
    
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result
                + ((mimatID == null) ? 0 : mimatID.hashCode());
        return result;
    }    
    
    
    /**
     * @return the mimatID
     */
    public String getMimatID() {
        return mimatID;
    }

    /**
     * @param mimatID the mimatID to set
     */
    public void setMimatID(String mimatID) {
        this.mimatID = mimatID;
    }

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * @param name the name to set
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @return the parent
     */
    public String getParent() {
        return parent;
    }

    /**
     * @param parent the parent to set
     */
    public void setParent(String parent) {
        this.parent = parent;
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
     * @return the sequence
     */
    public String getSequence() {
        return seq;
    }

    /**
     * @param sequence the sequence to set
     */
    public void setSequence(String sequence) {
        this.seq = sequence;
    }

    /**
     * @return the startPos
     */
    public int getStartPos() {
        return startPos;
    }

    /**
     * @param startPos the startPos to set
     */
    public void setStartPos(int startPos) {
        this.startPos = startPos;
    }

    /**
     * @return the endPos
     */
    public int getEndPos() {
        return endPos;
    }

    /**
     * @param endPos the endPos to set
     */
    public void setEndPos(int endPos) {
        this.endPos = endPos;
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

    /**
     * @return the isomiRString
     */
    public String getIsomiRString() {
        return isomiRString;
    }

    /**
     * @return the note
     */
    public String getNote() {
        return note;
    }

    /**
     * @param note the note to set
     */
    public void setNote(String note) {
        this.note = note;
    }

    /**
     * @return the host
     */
    public String getHost() {
        return host;
    }

    /**
     * @param host the host to set
     */
    public void setHost(String host) {
        this.host = host;
    }

    /**
     * @return the largestInternalLoop
     */
    public int getLargestInternalLoop() {
        return largestInternalBulge;
    }

    /**
     * @param largestInternalLoop the largestInternalLoop to set
     */
    public void setLargestInternalLoop(int largestInternalLoop) {
        this.largestInternalBulge = largestInternalLoop;
    }

    /**
     * @return the numOfInternalLoops
     */
    public int getNumOfInternalLoops() {
        return numOfBulges;
    }

    /**
     * @param numOfInternalLoops the numOfInternalLoops to set
     */
    public void setNumOfInternalLoops(int numOfInternalLoops) {
        this.numOfBulges = numOfInternalLoops;
    }

    /**
     * @return the numOfUnpairedBases
     */
    public int getNumOfUnpairedBases() {
        return numOfUnpairedBases;
    }

    /**
     * @param numOfUnpairedBases the numOfUnpairedBases to set
     */
    public void setNumOfUnpairedBases(int numOfUnpairedBases) {
        this.numOfUnpairedBases = numOfUnpairedBases;
    }

    /**
     * @return the fractOfUnpairedBases
     */
    public double getFractOfUnpairedBases() {
        return fractOfUnpairedBases;
    }

    /**
     * @param fractOfUnpairedBases the fractOfUnpairedBases to set
     */
    public void setFractOfUnpairedBases(double fractOfUnpairedBases) {
        this.fractOfUnpairedBases = fractOfUnpairedBases;
    }

    /**
     * @return the firstBase
     */
    public char getFirstBase() {
        return firstBase;
    }

    /**
     * @param firstBase the firstBase to set
     */
    public void setFirstBase(char firstBase) {
        this.firstBase = firstBase;
    }

    /**
     * @return the stability
     */
    public double getStability() {
        return stability;
    }

    /**
     * @param stability the stability to set
     */
    public void setStability(double stability) {
        this.stability = stability;
    }

    /**
     * @return the dangleBaseOne
     */
    public char getDangleBaseOne() {
        return dangleBaseOne;
    }

    /**
     * @param dangleBaseOne the dangleBaseOne to set
     */
    public void setDangleBaseOne(char dangleBaseOne) {
        this.dangleBaseOne = dangleBaseOne;
    }

    /**
     * @return the dangleBaseTwo
     */
    public char getDangleBaseTwo() {
        return dangleBaseTwo;
    }

    /**
     * @param dangleBaseTwo the dangleBaseTwo to set
     */
    public void setDangleBaseTwo(char dangleBaseTwo) {
        this.dangleBaseTwo = dangleBaseTwo;
    }

    /**
     * @return the miStart
     */
    public int getMiStart() {
        return miStart;
    }

    /**
     * @param miStart the miStart to set
     */
    public void setMiStart(int miStart) {
        this.miStart = miStart;
    }

    /**
     * @return the miEnd
     */
    public int getMiEnd() {
        return miEnd;
    }

    /**
     * @param miEnd the miEnd to set
     */
    public void setMiEnd(int miEnd) {
        this.miEnd = miEnd;
    }

    /**
     * @return the featureSet
     */
    public HashMap getFeatureSet() {
        return featureSet;
    }

    /**
     * @return the evidence
     */
    public String getEvidence() {
        return evidence;
    }

    /**
     * @param evidence the evidence to set
     */
    public void setEvidence(String evidence) {
        this.evidence = evidence;
    }

    /**
     * @return the references
     */
    public String getReferences() {
        return references;
    }

    /**
     * @param references the references to set
     */
    public void setReferences(String references) {
        this.references = references;
    }
    
}
