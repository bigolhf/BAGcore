/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core;

import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author weibo / simon rayner
 */
public class PreMiRNA extends SimpleRNASequence{

    public PreMiRNA(){
        super();
    }
    public PreMiRNA(String id,String seq){
        super(id,seq);
    }

    private String              host                            ="";
    private String              host3code                       ="";
    private String              note                           ="";
    private String              dbxrefs                         ="";
    
    private int                 maxInternalLoopSize             =0;
    private int                 numOfInternalLoops              =0;
    private int                 numOfUnpairedBases              =0;
    private double              fractOfUnpairedBases            =0;

    private int                 lowerStemSize                   =0;
    private int                 upperStemSize                   =0;
    private int                 topStemSize                     =0;
    private int                 upperStart                      =0;
    private int                 upperEnd                        =0;
    private int                 largestInternalLoopInLowerStem  =0;
    private int                 largestInternalLoopInTopStem    =0;
    private int                 numOfInternalLoopsInLowerStem   =0;
    private int                 numOfInternalLoopsInTopStem     =0;
    private int                 numOfUnpairedBasesInLowerStem   =0;
    private int                 numOfUnpairedBasesInTopStem     =0;
    private double              fractOfUnpairedBasesInLowerStem =0;
    private double              fractOfUnpairedBasesInTopStem   =0;

    private HashMap             featureSet;

    private ArrayList<MiRNAFeature>         miRNAList           = new ArrayList<>();
    private ArrayList<ShortPubMedEntry>     pubmedRefList       = new ArrayList<>();
    
    private int index=0;

    public void addProduct(MiRNAFeature miRNA){
        miRNAList.add(miRNA);
        index+=1;
    }

    /**
     * add pubmed reference to this pre-miRNA
     * 
     * @param newRef 
     */
    public void addPubmedRef(ShortPubMedEntry newRef){
        pubmedRefList.add(newRef);
    }

    /**
     * return number of miRNAs associated with this pre-miRNA
     * 
     * @return 
     */
    public int SizeOfProduct(){
        return miRNAList.size();
    }

    
    public MiRNAFeature NextMiRNA(){
        index-=1;
        return miRNAList.get(index);
    }

    
    /**
     * because the list of features we might calculate can change, we also store
     * in a hash map to make it simpler to collectively pass information
     */
    public void buildFeatureSet(){
        featureSet=new HashMap();
        
        featureSet.put("preRNA_sequence",               this.getSeq());
        featureSet.put("preRNA_structure",              this.getStructureStr());
        featureSet.put("preRNA_energy",                 this.getEnergy());
        featureSet.put("preRNA_size",                   this.getLength());
        featureSet.put("preRNA_GC_content",             this.getGC_content() );
        featureSet.put("preRNA_A_content",              this.getA_content());
        featureSet.put("preRNA_U_content",              this.getU_content());
        featureSet.put("preRNA_G_content",              this.getG_content());
        featureSet.put("preRNA_C_content",              this.getC_content());
        featureSet.put("preRNA_pair_number",            this.getNumberOfPairs());
        featureSet.put("preRNA_G-U_number",             this.getGU_num());
        featureSet.put("preRNA_unpair_number",          this.getUnpairedBase_num());
        featureSet.put("preRNA_unpair_rate",            this.getUnpairedBase_rate());
        featureSet.put("preRNA_internalLoop_number",    this.getInternalLoop_num());
        featureSet.put("preRNA_internalLoop_size",      this.getInternalLoopSize());
        featureSet.put("upperStem_start",               this.getUpperStart());
        featureSet.put("upperStem_end",                 this.getUpperEnd());
        featureSet.put("upperStem_size",                this.getUpperStemSize());
        featureSet.put("lowerStem_size",                this.getLowerStemSize());
        featureSet.put("lowerStem_unpair_number",       this.getLowerStemUnpairedBase_num());
        featureSet.put("lowerStem_unpair_rate",         this.getLowerStemUnpairedBase_rate());
        featureSet.put("lowerStem_internalLoop_number", this.getLowerStemInternalLoop_num());
        featureSet.put("lowerStem_internalLoop_size",   this.getLowerStemInternalLoopSize());
        featureSet.put("topStem_size",                  this.getTopStemSize());
        featureSet.put("topStem_unpair_number",         this.getTopStemUnpairedBase_num());
        featureSet.put("topStem_unpair_rate",           this.getTopStemUnpairedBase_rate());
        featureSet.put("topStem_internalLoop_number",   this.getTopStemInternalLoop_num());
        featureSet.put("topStem_internalLoop_size",     this.getTopStemInternalLoopSize());
    }

    public HashMap getFeatureSet(){
        return featureSet;
    }

    /**
     * @return the maxInternalLoopSize
     */
    public int getInternalLoopSize() {
        return maxInternalLoopSize;
    }

    /**
     * @param maxInternalLoopSize the maxInternalLoopSize to set
     */
    public void setBiggestInternalLoop(int maxInternalLoopSize) {
        this.maxInternalLoopSize = maxInternalLoopSize;
    }

    /**
     * @return the internalLoop_num
     */
    public int getInternalLoop_num() {
        return numOfInternalLoops;
    }

    /**
     * @param internalLoop_num the internalLoop_num to set
     */
    public void setNumberOfInternalLoops(int internalLoop_num) {
        this.numOfInternalLoops = internalLoop_num;
    }

    /**
     * @return the unpairedBase_num
     */
    public int getUnpairedBase_num() {
        return numOfUnpairedBases;
    }

    /**
     * @param unpairedBase_num the unpairedBase_num to set
     */
    public void setNumberOfUnpairedBases(int unpairedBase_num) {
        this.numOfUnpairedBases = unpairedBase_num;
    }

    /**
     * @return the unpairedBase_rate
     */
    public double getUnpairedBase_rate() {
        return fractOfUnpairedBases;
    }

    /**
     * @param unpairedBase_rate the unpairedBase_rate to set
     */
    public void setUnpairedBaseRate(double unpairedBase_rate) {
        this.fractOfUnpairedBases = unpairedBase_rate;
    }

     /**
     * @return the lowerStemSize
     */
    public int getLowerStemSize() {
        return lowerStemSize;
    }

    /**
     * @param lowerStemSize the lowerStemSize to set
     */
    public void setLowerStemLength(int lowerStemSize) {
        this.lowerStemSize = lowerStemSize;
    }

    /**
     * @return the upperStemSize
     */
    public int getUpperStemSize() {
        return upperStemSize;
    }

    /**
     * @param upperStemSize the upperStemSize to set
     */
    public void setUpperStemSize(int upperStemSize) {
        this.upperStemSize = upperStemSize;
    }

    /**
     * @return the topStemSize
     */
    public int getTopStemSize() {
        return topStemSize;
    }

    /**
     * @param topStemSize the topStemSize to set
     */
    public void setTopStemSize(int topStemSize) {
        this.topStemSize = topStemSize;
    }

    /**
     * @return the lowerStemInternalLoopSize
     */
    public int getLowerStemInternalLoopSize() {
        return largestInternalLoopInLowerStem;
    }

    /**
     * @param lowerStemInternalLoopSize the lowerStemInternalLoopSize to set
     */
    public void SetLargestInternalLoopInLowerStem(int lowerStemInternalLoopSize) {
        this.largestInternalLoopInLowerStem = lowerStemInternalLoopSize;
    }

    /**
     * @return the topStemInternalLoopSize
     */
    public int getTopStemInternalLoopSize() {
        return largestInternalLoopInTopStem;
    }

    /**
     * @param topStemInternalLoopSize the topStemInternalLoopSize to set
     */
    public void setTopStemInternalLoopSize(int topStemInternalLoopSize) {
        this.largestInternalLoopInTopStem = topStemInternalLoopSize;
    }

    /**
     * @return the lowerStemInternalLoop_num
     */
    public int getLowerStemInternalLoop_num() {
        return numOfInternalLoopsInLowerStem;
    }

    /**
     * @param lowerStemInternalLoop_num the lowerStemInternalLoop_num to set
     */
    public void getNumOfInternalLoopsInLowerStem(int lowerStemInternalLoop_num) {
        this.numOfInternalLoopsInLowerStem = lowerStemInternalLoop_num;
    }

    /**
     * @return the topStemInternalLoop_num
     */
    public int getTopStemInternalLoop_num() {
        return numOfInternalLoopsInTopStem;
    }

    /**
     * @param topStemInternalLoop_num the topStemInternalLoop_num to set
     */
    public void setNumOfInternalLoopsInTopStem(int topStemInternalLoop_num) {
        this.numOfInternalLoopsInTopStem = topStemInternalLoop_num;
    }

    /**
     * @return the lowerStemUnpairedBase_num
     */
    public int getLowerStemUnpairedBase_num() {
        return numOfUnpairedBasesInLowerStem;
    }

    /**
     * @param lowerStemUnpairedBase_num the lowerStemUnpairedBase_num to set
     */
    public void setNumOfUnpairedBasesInLowerStem(int lowerStemUnpairedBase_num) {
        this.numOfUnpairedBasesInLowerStem = lowerStemUnpairedBase_num;
    }

    /**
     * @return the topStemUnpairedBase_num
     */
    public int getTopStemUnpairedBase_num() {
        return numOfUnpairedBasesInTopStem;
    }

    /**
     * @param topStemUnpairedBase_num the topStemUnpairedBase_num to set
     */
    public void setNumOfUnpairedBasesInTopStem(int topStemUnpairedBase_num) {
        this.numOfUnpairedBasesInTopStem = topStemUnpairedBase_num;
    }

    /**
     * @return the lowerStemUnpairedBase_rate
     */
    public double getLowerStemUnpairedBase_rate() {
        return fractOfUnpairedBasesInLowerStem;
    }

    /**
     * @param lowerStemUnpairedBase_rate the lowerStemUnpairedBase_rate to set
     */
    public void setFractOfUnpairedBasesInLowerStem(double lowerStemUnpairedBase_rate) {
        this.fractOfUnpairedBasesInLowerStem = lowerStemUnpairedBase_rate;
    }

    /**
     * @return the topStemUnpairedBase_rate
     */
    public double getTopStemUnpairedBase_rate() {
        return fractOfUnpairedBasesInTopStem;
    }

    /**
     * @param topStemUnpairedBase_rate the topStemUnpairedBase_rate to set
     */
    public void setFractOfUnpairedBasesInTopStem(double topStemUnpairedBase_rate) {
        this.fractOfUnpairedBasesInTopStem = topStemUnpairedBase_rate;
    }

    public void setUpperStart(int upperStart) {
        this.upperStart = upperStart;
    }

    /**
     * @return the upperStart
     */
    public int getUpperStart() {
        return upperStart;
    }

    /**
     * @return the upperEnd
     */
    public int getUpperEnd() {
        return upperEnd;
    }

    /**
     * @param upperEnd the upperEnd to set
     */
    public void setUpperEnd(int upperEnd) {
        this.upperEnd = upperEnd;
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
     * @return the host3code
     */
    public String getHost3code() {
        return host3code;
    }

    /**
     * @param host3code the host3code to set
     */
    public void setHost3Lettercode(String host3code) {
        this.host3code = host3code;
    }

    /**
     * @return the notes
     */
    public String getNote() {
        return note;
    }

    /**
     * @param notes the notes to set
     */
    public void setNote(String notes) {
        this.note = notes;
    }

    /**
     * @return the dbxrefs
     */
    public String getDbxrefs() {
        return dbxrefs;
    }

    /**
     * @param dbxrefs the dbxrefs to set
     */
    public void setDbxrefs(String dbxrefs) {
        this.dbxrefs = dbxrefs;
    }



}

