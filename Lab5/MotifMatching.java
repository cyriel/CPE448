import java.util.*;
import java.io.*;

public class MotifMatching {
   static HashMap<Character, String> ref = new HashMap<Character, String>();
   static String fastafile, gfffile, motif, DNA = "";
   static ArrayList<Isoform> isoforms = new ArrayList<Isoform>();
   static ArrayList<DNAFrag> ups = new ArrayList<DNAFrag>();
   static ArrayList<Output> output = new ArrayList<Output>();
   static int match, miss, cutoff = 0;
   static int start = 0, end = -1, upstream;
   
   public static void hashInit() {
      ref.clear();
      ref.put('A',"A");
      ref.put('T',"T");
      ref.put('G',"G");
      ref.put('C',"C");
      ref.put('R',"AG");
      ref.put('Y',"CT");
      ref.put('W',"AT");
      ref.put('S',"CG");
      ref.put('M',"AC");
      ref.put('K',"GT");
      ref.put('H',"ACT");
      ref.put('B',"CGT");
      ref.put('V',"ACG");
      ref.put('D',"AGT");
      ref.put('N',"ACGT");
   }

   // Assuming we parsed dna to be the same length as motif
   // In the chance that it isn't, we will kill the program
   // I wish we could run fork bombs in languages like Java
   // Never say fork bomb on an airplane
   // Forks are very dangerous.
   public static int score(String dna) {
      int score = 0;
      for (int i = 0; i < dna.length(); i++)
         if (ref.get(motif.charAt(i)).indexOf(dna.charAt(i)) != -1)
            score += match;
         else
            score += miss;
      return score;
   }
   
   public static ArrayList<Output> search(DNAFrag frag) {
      int i, temp;
      ArrayList<Output> output = new ArrayList<Output>();
      
      for(i = frag.start; i+motif.length() < frag.stop; i = temp) {
         temp = i+motif.length();
         String seq = DNA.substring(i, temp);
         int value = score(seq);
         if (value > cutoff)
            output.add(new Output(frag.iso, seq, value, i, temp, start));
      }
      
      return output;
   }
   
   public static void purgeDuplicates() {
      for (int i = 0; i < output.size()-1; i++) {
         Output out1 = output.get(i);
         for (int n = i+1; n < output.size(); n++) {
            Output out2 = output.get(n);
            // if the start is the same for two outputs, then they have to be the same since the motif is the same
            if (out1.start == out2.start) {
               output.remove(n);
               n--;
            }
         }
      }
   }
   
   public static void parseConfigInput() {
      try {
         Scanner parser = new Scanner(new File("config.txt"));
         
         fastafile = parser.nextLine().split("=")[1].trim();
         gfffile = parser.nextLine().split("=")[1].trim();
         upstream = Integer.parseInt(parser.nextLine().split("=")[1].trim());
         motif = parser.nextLine().split("=")[1].trim();
         match = Integer.parseInt(parser.nextLine().split("=")[1].trim());
         miss = Integer.parseInt(parser.nextLine().split("=")[1].trim());
         cutoff = Integer.parseInt(parser.nextLine().split("=")[1].trim());
         
         /*
         System.out.println("FASTA file name = " + fastafile);
         System.out.println("GFF file name = " + gfffile);
         System.out.println("Upstream length = " + upstream);
         System.out.println("Motif = " + motif);
         System.out.println("Match = " + match);
         System.out.println("Miss = " + miss);
         System.out.println("Cutoff = " + cutoff);
         */
         
      } catch(FileNotFoundException e) {
         System.err.println("Config file not found!");
         System.exit(-1);
      } catch(NoSuchElementException ex) {
         System.err.println("Missing line in config file!");
         System.exit(-2);
      } catch(ArrayIndexOutOfBoundsException aioobe) {
         System.err.println("Missing equals sign in config file!");
         System.exit(-3);
      } catch(NumberFormatException nfe) {
         System.err.println("Invalid number format in one of the config file inputs!");
         System.exit(-4);
      }
   }
   
   public static void main(String[] args) {
      hashInit();
      parseConfigInput();
       
         try {
            Scanner fasta = new Scanner(new File(fastafile));
            fasta.nextLine();    // dump that junk line in the fasta file
            while (fasta.hasNextLine())
               DNA += fasta.nextLine();
//            System.out.println(DNA.length());
            if (end == -1)
               end = DNA.length();
            else
               DNA = DNA.substring(start,end);
            
            Scanner gff = new Scanner(new File(gfffile));
            
            Isoform isoform = new Isoform();
//            System.out.println("Done reading FASTA file");
            while (gff.hasNextLine()) {
               
               String tempGFF = gff.nextLine();
               Scanner parser = new Scanner(tempGFF);
               parser.next();
               parser.next();
               String feature = parser.next();
               int startRec = parser.nextInt()-1-start;
               int stopRec = parser.nextInt()-start;
               parser.next();
               boolean direction = parser.next().charAt(0) == '+' ? true: false;
               parser.next();
               parser.next();
               String isoName = parser.next();
               isoName = isoName.replace("\"","");
               isoName = isoName.replace(";","");
//               isoName = isoName.substring(0,isoName.indexOf('-'));
               if (!(startRec < 0 || stopRec > end - start)) {
                  if (isoform.name == "")    // initialization of isoform object
                     isoform.name = isoName;
                  else if (!isoform.name.substring(0,isoform.name.indexOf('-')).equals(isoName.substring(0,isoName.indexOf('-')))) { //add current isoform to ArrayList and create new isoform
                     if (isoform.start >= 0 && isoform.stop >=0) { // dump isoforms without a start and a stop
                        makeUps(isoform, upstream);
                        isoforms.add(isoform);
//                        System.out.println(isoform.name + " " + isoform.start + " " + isoform.stop);
                     }
                     isoform = new Isoform();
                     isoform.name = isoName;
                  }
                  // there is no else needed because the current isoform is still current (tautological?)
                  
                  if (isoform.name.equals(isoName)) {
                     Record record = new Record(feature, startRec, stopRec, direction, isoName);
      //               System.out.println(record.isoform + " " + record.feature);
                     
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
                     isoform.direction = direction;
                     isoform.records.add(record);
                  }
               }
            }
            if (!isoform.name.equals("")) {
               if (isoform.start >= 0 && isoform.stop >=0) {
//                  System.out.println(isoform.name + " " + isoform.start + " " + isoform.stop);
                  makeUps(isoform, upstream);
                  isoforms.add(isoform); // after we finish parsing GFF file, we should still have one isoform left to add to array list
               }
            }

            for(int i = 0; i<ups.size(); i++)
               output.addAll(ups.get(i).searched);
            
            purgeDuplicates();
            Collections.sort(output);
            
            File out = new File("Scores.csv");
            try {
               out.createNewFile();
               PrintWriter outPrinter = new PrintWriter(new FileWriter(out));
               
               outPrinter.println("Input Motif:," + motif);
               outPrinter.println("FASTA File Name:," + fastafile);
               outPrinter.println("Scoring Rules - match:," + match + ",miss:," + miss);
               outPrinter.println("Upstream Search Length:," + upstream);
               outPrinter.println("Cut-off Score:," + cutoff);
               outPrinter.println();
               outPrinter.println("Gene,Motif,Score,Start,Stop");
               
               for (int i = 0; i < output.size(); i++)
                  outPrinter.println(output.get(i));
               
               outPrinter.close();
            } catch (IOException e) {
               System.err.println("IO error occured");
            }

            
         } catch (FileNotFoundException e) {
            System.err.println("Input file not found");
         }
   }
   
   // need to change this to add to an ups array list
   public static void makeUps(Isoform iso, int size) {
      int start;
      DNAFrag frag;
      if (iso.direction) {
         start = iso.start - size;
         if (start < 0) start = 0;
         frag = new DNAFrag(start, iso.start);
         frag.iso = iso.name;
      } else {
         start = iso.stop + size;
         if (start > DNA.length()) start = DNA.length();
         frag = new DNAFrag(iso.stop, start);
         frag.iso = iso.name;
      }
//      System.out.println("frag start: " + frag.start + " stop: " + frag.stop);
      frag.searched = search(frag);
      if (!frag.searched.isEmpty())
         ups.add(frag);
   }
   
   // need a search function to search through the ups data and generate a tuple with gene upstream came from, sequence matched, score, and location of the sequence (start + stop)
}

class Isoform{
   ArrayList<Record> records = new ArrayList<Record>();
   String name = "";
   int start;
   int stop;
   boolean direction;
}

class Record{
   String feature, isoform;
   int start, stop;
   boolean direction; // true = forward, false = backwards
   
   public Record(String feature, int start, int stop, boolean direction, String isoform) {
      this.feature = feature;
      this.start = start;
      this.stop = stop;
      this.direction = direction;
      this.isoform = isoform;
   }
}

class DNAFrag {
   int start, stop;
   String iso;
   ArrayList<Output> searched;
   
   public DNAFrag(int start, int stop) {
      this.start = start;
      this.stop = stop;
   }
}

class Output implements Comparable<Output>{
   String gene, sequence;
   int score, start, stop;
   
   public Output(String gene, String sequence, int score, int start, int stop, int offset) {
      this.gene = gene;
      this.sequence = sequence;
      this.score = score;
      this.start = start+1+offset;
      this.stop = stop+offset;
   }
   
   public String toString() {
      return gene + "," + sequence + "," + score + "," + start + "," + stop; 
   }
   
   public int compareTo(Output other) {
      return (score == other.score) ? 0 : other.score - score;
   }
}
