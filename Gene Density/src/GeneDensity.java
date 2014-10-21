import java.util.*;
import java.io.*;

/*
 * 1 class per isoform (cut off everything after dash)
store list of isoform in arraylist

Class Name: isoform
Class Name: record
column 3: feature
column 4: start
column 5: stop
column 7: strand (boolean)
column 10: name (use function to clean "" and ""; && everything after dash)
 * 
 */

public class GeneDensity {
   static ArrayList<Isoform> isoforms = new ArrayList<Isoform>(); //place isoform objects here
   static ArrayList<DNAFrag> nests = new ArrayList<DNAFrag>(); // place discovered nests here
   static int length, numA, numT, numG, numC, numGC, numN;
   static int nestCount=0, nestLength=0, nestA=0, nestT=0, nestG=0, nestC=0, nestGC=0, nestN=0;
   static int igCount=0, igLength=0, igA, igT, igG, igC, igGC, igN;
   static int upsCount=0, upsLength=0, upsA=0, upsT=0, upsG=0, upsC=0, upsGC=0, upsN=0;
   static boolean isoformSorted = false;
   static String DNA;
   static int[] exon = new int[8], cds = new int[8], intron = new int[8], gene = new int[8];
   static int start = 1, upstream = 0, end = -1;    // default start, stop, and upstream values
   
   // need to implement a startOffset and an endOffset to handle searching through the GFF file for applicable ranges
   // this /could/ still allow upstream to be calculated at full distances
   public static void main(String[] args) {
      String config = "config.txt";
      Scanner parser;
      try {
         parser = new Scanner(new File(config));
         String fastafile = parser.nextLine();
         String gfffile = parser.nextLine();
         if (parser.hasNextInt()) {
            upstream = parser.nextInt();
            if (parser.hasNextInt()) {
               start = parser.nextInt();
               if (parser.hasNextInt())
                  end = parser.nextInt();
               else
                  start = 1;
            }
         }
         start--;
         try {
            Scanner fasta = new Scanner(new File(fastafile));
            fasta.nextLine();    // dump that junk line in the fasta file
            DNA = "";     // initialize dna string
            while (fasta.hasNextLine())
               DNA += fasta.nextLine();
            if (end == -1)
               end = DNA.length();
            String dna_fragment = DNA.substring(start,end);
            
            Scanner gff = new Scanner(new File(gfffile));
            
            Isoform isoform = new Isoform();
            while (gff.hasNextLine()) {
               String tempGFF = gff.nextLine();
               parser = new Scanner(tempGFF);
               parser.next();
               parser.next();
               String feature = parser.next();
               int startRec = parser.nextInt()-1;
               int stopRec = parser.nextInt();
               parser.next();
               boolean direction = parser.next().charAt(0) == '+' ? true: false;
               parser.next();
               parser.next();
               String isoName = parser.next();
               isoName = isoName.replace("\"","");
               isoName = isoName.substring(0,isoName.indexOf('-'));
               
               if (isoform.name == "")    // initialization of isoform object
                  isoform.name = isoName;
               else if (!isoform.name.equals(isoName)) { //add current isoform to ArrayList and create new isoform
                  isoform.runAllDataCalc(extractGeneDNA(DNA,isoform));
                  calcUpsData(isoform, upstream);
                  isoforms.add(isoform);
//                  System.out.println(isoform.name + " " + isoform.start + " " + isoform.stop + " " + isoform.length);
                  isoform = new Isoform();
                  isoform.name = isoName;
               }
               // there is no else needed because the current isoform is still current (tautological?)
               
               Record record = new Record(feature, startRec, stopRec, direction, isoName);
//               System.out.println(record.isoform + " " + record.feature);
               record.genData(extractRecordDNA(DNA,record));
               
               if (record.feature.equals("start_codon"))
                  if (direction)
                     isoform.start = startRec;
                  else
                     isoform.stop = stopRec;
               else if (record.feature.equals("stop_codon"))
                  if (direction)
                     isoform.stop = stopRec;
                  else
                     isoform.start = startRec;
               else if (record.feature.equals("exon"))       // not really necessary, but don't want to chance a random feature name from getting input
                  isoform.exons.add(record);
               isoform.direction = direction;
               isoform.records.add(record);
            }
            isoform.runAllDataCalc(extractGeneDNA(DNA,isoform));
//            System.out.println(isoform.name + " " + isoform.start + " " + isoform.stop + " " + isoform.length);
            calcUpsData(isoform, upstream);
            isoforms.add(isoform); // after we finish parsing GFF file, we should still have one isoform left to add to array list
            
            genData(DNA);
            calcNests();
            calcIntergenic();
            
            genCSVdata();
            genCSVnests();
            
         } catch (FileNotFoundException e) {
            System.err.println("Input file not found");
         }
      } catch (FileNotFoundException e) {
         System.err.println("Config file not found");
      }
      
   }

   public static String toCSV(int[] input) {
      String returned = "" + input[0];
      for (int i = 1; i < input.length; i++)
         returned += "," + input[i];
      return returned;
   }
   
   public static int[] updateArray(int[] array, Record rec) {
      int[] returned = array;
      returned[0]++;
      returned[1] += rec.length;
      returned[2] += rec.numA;
      returned[3] += rec.numT;
      returned[4] += rec.numG;
      returned[5] += rec.numC;
      returned[6] += rec.numGC;
      returned[7] += rec.numN;
      return returned;
   }
   
   // simultaneously generate CSV data for genes and the records in the gene
   public static void genCSVdata() {
      File gen = new File("GeneralData.csv");
      File iso = new File("GeneData.csv");
      double[] geneData = new double[7];
      try {
         gen.createNewFile();
         iso.createNewFile();
         PrintWriter printerGen = new PrintWriter(new FileWriter(gen));
         PrintWriter printerIso = new PrintWriter(new FileWriter(iso));
         printerGen.println("Feature,Count,Length,%A,%T,%G,%C,%GC,%N");
         printerIso.println("Gene,CDS Count,CDS Length,Exon Count,Exon Length,Intron Count,Intron Length,Gene Length");
         for (Isoform isoform: isoforms) {
            for (Record rec: isoform.records) {
               if (rec.feature.equals("exon"))
                  exon = updateArray(exon, rec);
               else if (rec.feature.equals("CDS"))
                  cds = updateArray(cds, rec);
            }
            for (int i = 0; i < 8; i++)
               intron[i] += isoform.intron[i];
            gene[0]++;
            gene[1] += isoform.length;
            gene[2] += isoform.numA;
            gene[3] += isoform.numT;
            gene[4] += isoform.numG;
            gene[5] += isoform.numC;
            gene[6] += isoform.numGC;
            gene[7] += isoform.numN;
            
            geneData[0] += isoform.CDScount;
            geneData[1] += isoform.CDSlength;
            geneData[2] += isoform.exonCount;
            geneData[3] += isoform.exonLength;
            geneData[4] += isoform.intronCount;
            geneData[5] += isoform.intronLength;
            geneData[6] += isoform.length;
            
            printerIso.println(isoform.name+","+isoform.CDScount+","+isoform.CDSlength+","+isoform.exonCount+","+isoform.exonLength+","+isoform.intronCount+","+isoform.intronLength+","+isoform.length);
         }
         printerGen.println("Exons,"+exon[0]+","+exon[1]+","+(double)exon[2]/exon[1]+","+(double)exon[3]/exon[1]+","+(double)exon[4]/exon[1]+","+(double)exon[5]/exon[1]+","+(double)exon[6]/exon[1]+","+(double)exon[7]/exon[1]);
         printerGen.println("CDS,"+cds[0]+","+cds[1]+","+(double)cds[2]/cds[1]+","+(double)cds[3]/cds[1]+","+(double)cds[4]/cds[1]+","+(double)cds[5]/cds[1]+","+(double)cds[6]/cds[1]+","+(double)cds[7]/cds[1]);
         printerGen.println("Intron,"+intron[0]+","+intron[1]+","+(double)intron[2]/intron[1]+","+(double)intron[3]/intron[1]+","+(double)intron[4]/intron[1]+","+(double)intron[5]/intron[1]+","+(double)intron[6]/intron[1]+","+(double)intron[7]/intron[1]);
         printerGen.println("Intergenic Regions,"+igCount+","+igLength+","+(double)igA/igLength+","+(double)igT/igLength+","+(double)igG/igLength+","+(double)igC/igLength+","+(double)igGC/igLength+","+(double)igN/igLength);
         printerGen.println("Genes,"+gene[0]+","+gene[1]+","+(double)gene[2]/gene[1]+","+(double)gene[3]/gene[1]+","+(double)gene[4]/gene[1]+","+(double)gene[5]/gene[1]+","+(double)gene[6]/gene[1]+","+(double)gene[7]/gene[1]);
         printerGen.println("Upstream (" + upstream + " b),"+upsCount+","+upsLength+","+(double)upsA/upsLength+","+(double)upsT/upsLength+","+(double)upsG/upsLength+","+(double)upsC/upsLength+","+(double)upsGC/upsLength+","+(double)upsN/upsLength);
         printerGen.println("Nested Genes,"+nestCount+","+nestLength+","+(double)nestA/nestLength+","+(double)nestT/nestLength+","+(double)nestG/nestLength+","+(double)nestC/nestLength+","+(double)nestGC/nestLength+","+(double)nestN/nestLength);
         printerGen.close();
         for (int i = 0; i < 7; i++)
            geneData[i] = 1.0*geneData[i]/isoforms.size();
         printerIso.println("AVERAGE"+","+geneData[0]+","+geneData[1]+","+geneData[2]+","+geneData[3]+","+geneData[4]+","+geneData[5]+","+geneData[6]);
         printerIso.close();
      } catch (IOException e) {
         System.err.println("IO error occured");
      }

   }
   
   public static void genCSVnests() {
      File out = new File("nests.csv");
      try {
         out.createNewFile();
         FileWriter writer = new FileWriter(out);
         PrintWriter printer = new PrintWriter(writer);
         if (nests.size() == 0)
            printer.println("No nested genes were found");
         else {
            printer.println("Isoform 1,Isoform 2,Start,Stop");
            for (int i = 0; i < nests.size(); i++)
               printer.println(nests.get(i).toCSV());
         }
         printer.close();
      } catch (IOException e) {
         System.err.println("IO error occured");
      }
   }
   
   public static void calcUpsData(Isoform iso, int size) {
      int start;
      DNAFrag frag;
      if (iso.direction) {
         start = iso.start - size;
         if (start < 0) start = 0;
         frag = new DNAFrag(start, iso.start);
      } else {
         start = iso.stop + size;
         if (start > DNA.length()) start = DNA.length();
         frag = new DNAFrag(iso.stop, start);
      }
      frag.genData(extractFragDNA(DNA, frag));
      if (frag.length > 0) {
         upsCount++;
         upsLength += frag.length;
         upsA += frag.numA;
         upsT += frag.numT;
         upsG += frag.numG;
         upsC += frag.numC;
         upsGC += frag.numGC;
         upsN += frag.numN;
      }
   }
   
   // put this into Record's genData function to get statistical calculations
   public static String extractRecordDNA(String dna, Record rec) {
      String temp = dna.substring(rec.start, rec.stop);
      return rec.direction ? temp : reverseComplement(temp);
   }
   
   // this is here just because of the other ones. gotta keep things looking organized
   public static String extractFragDNA(String dna, DNAFrag frag) {
      return dna.substring(frag.start, frag.stop);
   }
   
   // not sure what this should be used for yet
   public static String extractGeneDNA(String dna, Isoform iso) {
      String temp = dna.substring(iso.start, iso.stop);
//      System.out.println(iso.start + " " + iso.stop + " " + temp.length() + " "  + reverseComplement(temp).length());
      return iso.direction ? temp : reverseComplement(temp);
   }
   
   // gets the reverse complement of a dna string. useful when the gene is written reversed
   public static String reverseComplement(String dna) {
      char[] reverse = new StringBuilder(dna).reverse().toString().toCharArray();
      String revCom = "";
      for (int i = 0; i<reverse.length; i++) {
         if (reverse[i] == 'A')
            revCom += 'T';
         else if (reverse[i] == 'T')
            revCom += 'A';
         else if (reverse[i] == 'G')
            revCom += 'C';
         else if (reverse[i] == 'C')
            revCom += 'G';
      }
      return revCom;
   }
   
   /* return values: 0: no nesting
    *                1: nested from iso1.start to iso2.stop
    *                2: nested from iso2.start to iso1.start
    *                3: iso1 is nested inside iso2
    *                4: nested from iso1.stop to iso2.stop
    *                5: nested from iso2.start to iso1.stop
    *                6: iso1 is nested inside iso2
    *                7: nested from iso2.start to iso1.stop
    *                8: nested from iso1.start to iso2.start
    *                9: iso2 is nested inside iso1
    *               10: nested from iso2.stop to iso1.stop
    *               11: nested from iso1.start to iso2.stop
    *               12: iso2 is nested inside iso1
    */
   public static int isNested(Isoform iso1, Isoform iso2) {
      if (iso1.start > iso2.start && iso1.start < iso2.stop)
         if (iso1.stop > iso2.stop)                   // nested from iso1.start to iso2.stop
            return 1;
         else if (iso1.stop < iso2.start)             // nested from iso2.start to iso1.start
            return 2;
         else                                         // iso1 is nested inside iso2 <-- something went wrong?
            return 3;
      else if (iso1.stop > iso2.start && iso1.stop < iso2.stop)
         if (iso1.start > iso2.stop)                  // nested from iso1.stop to iso2.stop
            return 4;
         else if (iso1.start < iso2.start)            // nested from iso2.start to iso1.stop
            return 5;
         else                                         // iso1 is nested inside iso2 <-- something went wrong?
            return 6;
      else if (iso2.start > iso1.start && iso2.start < iso1.stop)
         if (iso2.stop > iso1.stop)                   // nested from iso2.start to iso1.stop
            return 7;
         else if (iso2.stop < iso1.start)             // nested from iso1.start to iso2.start
            return 8;
         else                                         // iso2 is nested inside iso1 <-- something went wrong?
            return 9;
      else if (iso2.stop > iso1.start && iso2.stop < iso1.stop)
         if (iso2.start > iso1.stop)                  // nested from iso2.stop to iso1.stop
            return 10;
         else if (iso2.start < iso1.start)            // nested from iso1.start to iso2.stop
            return 11;
         else                                         // iso2 is nested inside iso1 <-- something went wrong?
            return 12;
      else
         return 0;
   }
   
   // hint should be return value from isNested
   // return is [start,stop]
   public static int[] nested(Isoform iso1, Isoform iso2, int hint) {
      int[] range;
      if (hint == 1)
         range = new int[] {iso1.start,iso2.stop};
      else if (hint == 2)
         range = new int[] {iso2.start,iso1.start};
      else if (hint == 4)
         range = new int[] {iso1.stop,iso2.stop};
      else if (hint == 5)
         range = new int[] {iso2.start,iso1.stop};
      else if (hint == 3 || hint == 6)
         range = new int[] {iso1.start,iso1.stop};
      else if (hint == 7)
         range = new int[] {iso2.start,iso1.stop};
      else if (hint == 8)
         range = new int[] {iso1.start,iso2.start};
      else if (hint == 10)
         range = new int[] {iso2.stop,iso1.stop};
      else if (hint == 11)
         range = new int[] {iso1.start,iso2.stop};
      else if (hint == 9 || hint == 12)
         range = new int[] {iso2.start,iso2.stop};
      else
         range = new int[] {0,0};
         
      return range;
   }
   
   // do an n^2 comparison of isoforms and find nested regions and store it as a DNAFrag in nests
   // finds and calculates Nests/data
   public static void calcNests() {
      // for each nested region, create a DNAFrag object, calculate the data, and then store it in the array list
      for (int i = 0; i < isoforms.size()-1; i++)
         for (int n = 1; n < isoforms.size(); n++) {
            int hint = isNested(isoforms.get(i), isoforms.get(n));
            if (hint != 0) {
               int[] startstop = nested(isoforms.get(i), isoforms.get(n), hint);
               DNAFrag frag = new DNAFrag(startstop[0], startstop[1]);
//               System.out.println(isoforms.get(i).name + "\t" + isoforms.get(n).name);
//               System.out.println(hint + "\t" + startstop[0] + "\t" + startstop[1]);
               frag.genData(extractFragDNA(DNA, frag));
               frag.iso1 = isoforms.get(i).name;
               frag.iso2 = isoforms.get(n).name;
               nestLength += frag.length;
               nestA += frag.numA;
               nestT += frag.numT;
               nestG += frag.numG;
               nestC += frag.numC;
               nestGC += frag.numGC;
               nestN += frag.numN;
               nests.add(frag);
            }
         }
      nestCount=nests.size();
   }
   
   // we can't generate nested and intergenic at the same time in the /veryveryvery/ rare case where we have multiple nests involving the same few genes
   public static void calcIntergenic() {
      // need to sort Isoforms from lowest to highest and look at the gaps between neighbors.
      // this function is going to pretend that overlaps don't happen, so stop is between another's start and stop, we ignore it.
      // ignore the previous line....if there is an overlap, we would be able to fix that by adding in the nested regions
      if (!isoformSorted)
         sortIsoforms();
      
      igLength = length;
      igA = numA;
      igT = numT;
      igG = numG;
      igC = numC;
      igGC = numGC;
      igN = numN;
      igCount = isoforms.size()+1; // there is a good chance that this logic is correct. if we really want to be accurate and test all cases, we should have two booleans in the main to verify start and stop location

      for(Isoform iso: isoforms) {
         igLength -= iso.length;
         igA -= iso.numA;
         igT -= iso.numT;
         igG -= iso.numG;
         igC -= iso.numC;
         igGC -= iso.numGC;
         igN -= iso.numN;
      }
      
   }
   
   public static void sortIsoforms() {
      Collections.sort(isoforms);
      isoformSorted = true;
   }
   
   public static void genData(String dna) {
      length = dna.length();
      numA = length - dna.replace("A", "").length();
      numT = length - dna.replace("T", "").length();
      numG = length - dna.replace("G", "").length();
      numC = length - dna.replace("C", "").length();
      numN = length - dna.replace("N", "").length();
      numGC = numG + numC;
   }
   
}

class Isoform implements Comparable<Isoform>{
   ArrayList<Record> records = new ArrayList<Record>();
   ArrayList<Record> exons = new ArrayList<Record>();
   String name = "";
   int start;
   int stop;
   boolean direction, exonSorted = false;
   int length, numA, numT, numG, numC, numGC, numN;
   int[] intron = new int[8];
   int CDScount=0, CDSlength=0, exonCount=0, exonLength=0, intronCount, intronLength;
   
   public void runAllDataCalc(String dna) {
      genData(dna);
      calcIntronData();
      genRecordStat();
   }
   
   public int compareTo(Isoform other) {
      if (start == other.start)
         return (stop == other.stop) ? 0 : stop - other.stop;
      else
         return start - other.start;
   }
   
   public void calcIntronData() {
      if (!exonSorted)
         sortExons();
      
      intron[0] = intronCount = exons.size()-1;
      intron[1] = length;
      intron[2] = numA;
      intron[3] = numT;
      intron[4] = numG;
      intron[5] = numC;
      intron[6] = numGC;
      intron[7] = numN;

      for(Record exon: records) {
         intron[1] -= exon.length;
         intron[2] -= exon.numA;
         intron[3] -= exon.numT;
         intron[4] -= exon.numG;
         intron[5] -= exon.numC;
         intron[6] -= exon.numGC;
         intron[7] -= exon.numN;
      }
      intronLength = intron[1];
   }
   
   // data pertaining to the gene as a whole, not the Records
   // assume that the dna string is pre-parsed using extractIsoformDNA()
   public void genData(String dna) {
      length = dna.length();
      numA = length - dna.replace("A", "").length();
      numT = length - dna.replace("T", "").length();
      numG = length - dna.replace("G", "").length();
      numC = length - dna.replace("C", "").length();
      numN = length - dna.replace("N", "").length();
      numGC = numG + numC;
   }
   
   public void genRecordStat() {
      for (Record rec: records) {
         if (rec.type == 3) {
            CDScount++;
            CDSlength += rec.length;
         } else if (rec.type == 4) {
            exonCount++;
            exonLength += rec.length;
         }
      }
   }
   
   public void sortExons() {
      Collections.sort(records);
      
      exonSorted = true;
   }
}

class Record implements Comparable<Record>{
   String feature, isoform;
   int start, stop;
   boolean direction; // true = forward, false = backwards
   int type = 0; //0: no type, 1: start, 2: stop, 3: CDS, 4: exon
   int length, numA, numT, numG, numC, numGC, numN;
   
   public Record(String feature, int start, int stop, boolean direction, String isoform) {
      this.feature = feature;
      this.start = start;
      this.stop = stop;
      this.direction = direction;
      this.isoform = isoform;
      
      if (feature.equals("start_codon"))
         type = 1;
      else if (feature.equals("stop_codon"))
         type = 2;
      else if (feature.equals("CDS"))
         type = 3;
      else if (feature.equals("exon"))
         type = 4;
   }
   
   public int compareTo(Record other) {
      if (start == other.start)
         return (stop == other.stop) ? 0 : stop - other.stop;
      else
         return start - other.start;
   }
   
   // assume that the dna string is pre-parsed using extractRecordDNA()
   // hopefully this is run before it becomes a part of an Isoform otherwise running is a pain
   public void genData(String dna) {
      length = dna.length();
      numA = length - dna.replace("A", "").length();
      numT = length - dna.replace("T", "").length();
      numG = length - dna.replace("G", "").length();
      numC = length - dna.replace("C", "").length();
      numN = length - dna.replace("N", "").length();
      numGC = numG + numC;
   }
}

class DNAFrag {
   int start, stop;
   int length, numA, numT, numG, numC, numGC, numN;
   String iso1, iso2; // this is where we put either the two isoforms that are nested or the isoforms surrounding the intergenic region
   
   public DNAFrag(int start, int stop) {
      this.start = start;
      this.stop = stop;
   }
   
   public String toCSV() {
      return iso1 + "," + iso2 + "," + start + "," + stop;
   }
   
   // assume that the dna string is pre-parsed using extractFragDNA()
   public void genData(String dna) {
      length = dna.length();
      numA = length - dna.replace("A", "").length();
      numT = length - dna.replace("T", "").length();
      numG = length - dna.replace("G", "").length();
      numC = length - dna.replace("C", "").length();
      numN = length - dna.replace("N", "").length();
      numGC = numG + numC;
   }
}
