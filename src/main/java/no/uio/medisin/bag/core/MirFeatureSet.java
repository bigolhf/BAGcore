/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core;

import no.uio.medisin.bag.core.MiRNAFeature;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 * Stores a complete miRBase release of miRNA data
 * 
 * @author sr
 */
public class MirFeatureSet {
    
    static  Logger                      logger                              = LogManager.getLogger();
    
    private List<MiRNAFeature>          miRBaseMiRNAList                    = new ArrayList<>();

    
    
    /**
     * 1. 
     * load miRNA specs (name and Chromosome position) from GFF file
     * downloaded from miRBase
     * Because different releases of miRBase use different releases of
     * reference genome, we have to track both miRBase and genome reference IDs
     * 2.
     * Load Sequence from Mature.fa file
     * 
     * @param host
     * @param gffFilePath           
     * @param fastaFilePath         
     * @throws IOException
     * 
     */
    public void loadMiRBaseData(String host, String gffFilePath, String fastaFilePath) throws IOException{

        HashMap <String, String> miRBaseSeq = new HashMap();
        try{
            logger.info("reading miRBase fasta reference file <" +  fastaFilePath + ">");
            BufferedReader brFA = new BufferedReader(new FileReader(new File(fastaFilePath)));
            String lineFA = null;
            while ((lineFA = brFA.readLine())!=null){

                String seq = brFA.readLine().trim();
                String entryHost = lineFA.split(" ")[0].substring(1).split("-")[0].trim();
                if(entryHost.equals(host)){
                    String mimatID = lineFA.split(" ")[1].trim();
                    miRBaseSeq.put(mimatID, seq);
                }

            }
            brFA.close();
            logger.info("read " + miRBaseSeq.size() + " entries");
            logger.info("done");

        }
        catch(IOException ex){
            logger.error("error reading miRBase fasta reference file <" +  fastaFilePath + ">\n" + ex.toString());
            throw new IOException("error reading miRBase fasta reference file <" +  fastaFilePath + ">");
        }
        
        
        
        try{
            String line = null;
            logger.info("reading miRBase fasta reference file <" +  gffFilePath + ">");
            BufferedReader brMiR = new BufferedReader(new FileReader(new File(gffFilePath)));
                while((line = brMiR.readLine())!= null){

                    if(line.startsWith("#")) continue;
                    if(line.contains("miRNA_primary_transcript")) continue;
                    /*
                        chr1            chromosome
                        source          n/a here
                        miRNA           feature type (n/a)
                        start pos
                        end pos
                        score           n/a here               
                        strand          (+/-)
                        frame           n/a here
                        attributes      e.g. ID=MIMAT0027619;Alias=MIMAT0027619;Name=hsa-miR-6859-3p;Derives_from=MI0022705

                    */
                    String chr = line.split("\t")[0].trim();
                    if(chr.contains("chr")) chr = chr.replace("chr", "");
                    int startPos = Integer.parseInt(line.split("\t")[3].trim());
                    int endPos = Integer.parseInt(line.split("\t")[4].trim());
                    String strand = line.split("\t")[6].trim();

                    String id = "";
                    String name = "";
                    String parent = "";
                    String attribs[] = line.split("\t")[8].split(";");

                    for (String attribStr: attribs){
                        String attribType = attribStr.split("=")[0].trim();
                        String attribValue = attribStr.split("=")[1].trim();
                        switch (attribType){
                            case "ID":
                                id = attribValue;
                                break;

                            case "Alias":                            
                                break;

                            case "Name":
                                name = attribValue;
                                break;

                            case "Derives_from":
                                parent = attribValue;
                                break;

                            default:
                                logger.warn("unknown attribute in parsing miRNA entry from GFF file " + gffFilePath + ">");
                                break;
                        }
                    }
                    // this may need revising. In some cases we dont care if there is no sequence, we still want the feature
                    String seq = miRBaseSeq.get(id);
                    this.miRBaseMiRNAList.add(new MiRNAFeature(name, chr, startPos, endPos, strand, id, parent, seq));
                    if(seq == null) 
                        logger.warn("no sequence found for entry <" + id + ">. Skipping");

                }
            brMiR.close();
            logger.info("read " + miRBaseMiRNAList.size() + " miRNA entries");
        }
        catch(IOException ex){
            logger.error("error reading miRBase gff reference file <" +  gffFilePath + ">\n" + ex.toString());
            throw new IOException("error reading miRBase gff reference file <" +  gffFilePath + ">");
        }
        
    }
    

    
    
    
    /**
     * Does the read sufficiently overlap a defined miRNA entry?
     * 
     * @param start
     * @param stop
     * @param chr
     * @param strand
     * @param bleed     int : specifies how much a read can 'miss' an entry
     *                        and still be counted
     * 
     * @return MiRNAFeature
     */
    public MiRNAFeature doesReadOverlapKnownMiRNA(int start, int stop, String chr, String strand, int bleed){
        
        for(MiRNAFeature miRBaseEntry: this.getMiRBaseMiRNAList()){
            if (miRBaseEntry.chromosomeMatch(chr)){
                if(strand.equals(miRBaseEntry.getStrand())){
                    
                    if( java.lang.Math.abs(start - miRBaseEntry.getStartPos()) <= bleed){

                        if( java.lang.Math.abs(stop - miRBaseEntry.getEndPos()) <= bleed){
                            return miRBaseEntry;
                        }

                    }

                }
            }
        }
        
        return null;
        
    }
    
    
    
    
    
    
    /**
     * get the number of entries in this release
     * 
     * @return 
     */
    public int getNumberOfEntries(){
        return getMiRBaseMiRNAList().size();
    }

    /**
     * @return the miRBaseMiRNAList
     */
    public List<MiRNAFeature> getMiRBaseMiRNAList() {
        return miRBaseMiRNAList;
    }
    
}
