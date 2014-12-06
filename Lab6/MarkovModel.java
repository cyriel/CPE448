import java.util.*;
import java.io.*;

public class MarkovModel {
   static String fastaFile;
//   static String gffFile;
   static int start;
   static int end;
   static float[][] transition;
   static float[][] emission;
   static String frag;
   
   public static final int SS31 = 0;
   public static final int SS32 = 1;
   public static final int I = 2;
   public static final int SS51 = 3;
   public static final int SS52 = 4;
   public static final int E = 5;
   public static final int NUMSTATES = 6;
   
   public static final int A = 0;
   public static final int T = 1;
   public static final int C = 2;
   public static final int G = 3;
   
   public static int size;
   
   public static String[] probState;
   public static int[] pathState;
   public static int probStart, probStop;
   public static float probPath;
   
   // observations needs to be the nucleotide sequence
   //public static final int OBSERVATIONS = {A, T, C, G};
   public static final int[] STATES = {SS31, SS32, I, SS51, SS52, E};
   
   public static float[][] viterbi;
   
   public static void main(String[] args) {
      transition = new float[6][6];
      emission = new float[6][4];
      parseConfigInput();
      // need to read the dna to variable. Look at previous labs for this
      size = end - start;
      viterbi = new float[size][NUMSTATES];
      probState = new String[size];
      pathState = new int[size];
      
      try {
         Scanner fasta = new Scanner(new File(fastaFile));
         fasta.nextLine();
         String dna = "";
         while (fasta.hasNextLine())
            dna += fasta.nextLine();
         frag = dna.substring(start,end);
         // this is where the magic needs to happen
         
         viterbiFunct(transition, emission);
         
         print();
         
         File txtOut = new File("output.txt");
         try {
            txtOut.createNewFile();
            PrintWriter txtPrinter = new PrintWriter(new FileWriter(txtOut));
            //String.format("%.2f", 100*probA)
            probStop = -1; // if the stop was not detected, report as error
            probStart = -1; // if the start was not detected, report error
            for (int i = 0; i < size; i++)
               if (i > 0 && probStart == -1 && pathState[i-1] == SS52 && pathState[i] == E)
                  probStart = i;
               else if (i < size-1 && probStop == -1 && pathState[i+1] == SS31 && pathState[i] == E)
                  probStop = i;
            if (probStart == -1)
               txtPrinter.println("Exon was not found.");
            else {
               txtPrinter.println("Start: " + (start + probStart+1));
               if (probStop == -1)
                  txtPrinter.println("Exon stop was not found.");
               else
                  txtPrinter.println("Stop: " + (start + probStop+1));
            }
            txtPrinter.println("Probability: " + String.format("%.4f", probPath));
            txtPrinter.println();
            txtPrinter.println("Most probable path:");
            for (int i = 0; i < size; i++)
               txtPrinter.println("\t" + stateVal(pathState[i]));
            txtPrinter.close();
         } catch (IOException e) {
            System.err.println("IO error occured");
         }
      } catch (FileNotFoundException e) {
         System.err.println("FASTA file not found");
      }
   }
   
   public static void print() {
      for (int n = 0; n < NUMSTATES; n++)
         System.out.print(stateVal(n) + "\t\t");
      System.out.println();
      for (int i = 0; i < 10; i++) {
         for (int n = 0; n < NUMSTATES; n++)
            System.out.print(viterbi[i][n] + "\t\t");
         System.out.println();
      }
   }
   
   public static void parseConfigInput() {
      try {
         Scanner parser = new Scanner(new File("config.txt"));
      
         fastaFile = parser.nextLine().split("=")[1].trim();
//         gffFile = parser.nextLine().split("=")[1].trim();
         start = Integer.parseInt(parser.nextLine().split("=")[1].trim());
         start--;
         end = Integer.parseInt(parser.nextLine().split("=")[1].trim());

         parser.nextLine(); parser.nextLine(); parser.nextLine(); //skip blank line and table labels
      
         // parse transition probability table
         for(int i = 0; i < 6; i++) {
            parser.next(); // discard row label
            for(int j = 0; j < 6; j++) {
               transition[i][j] = parser.nextFloat();
               //System.out.println("transition[" + i + "][" + j + "] = " + transition[i][j]);
            }
         }

         parser.nextLine(); parser.nextLine(); parser.nextLine(); parser.nextLine();//skip blank line and table label

         // parse emission probability table
         for(int i = 0; i < 6; i++) {
            parser.next(); // discard row label
            for(int j = 0; j < 4; j++) {
               emission[i][j] = parser.nextFloat();
               //System.out.println("emission[" + i + "][" + j + "] = " + emission[i][j]);
            }
         }
         /*
         System.out.println("FASTA file name = " + fastaFile);
         System.out.println("GFF file name = " + gffFile);
         System.out.println("Start coordinate = " + start);
         System.out.println("End coordinate = " + end);
         System.out.println("Transition[I][5'SS1] = " + transition[I][SS51]);
         System.out.println("Emission[SS52][C] = " + emission[SS52][C]);
         */
      
      } catch(FileNotFoundException e) {
         System.out.println("Config file not found!");
         return;
      } catch(NoSuchElementException ex) {
         System.out.println("Missing line in config file!");
         return;
      } catch(ArrayIndexOutOfBoundsException aioobe) {
         System.out.println("Missing equals sign in config file!");
         return;
      } catch(NumberFormatException nfe) {
         System.out.println("Invalid number format in one of the config file inputs!");
         return;
      }
   }
   
   public static void viterbiFunct(float[][] t_prob, float[][] e_prob) {
      viterbi[0][0] = 1; // everything is initialized as 0 so only need to set this for initialization
      
      for (int t = 1; t < size; t++) {
         float max;
         // P[] = e_prob
         // a[] = t_prob
         // n = NUMSTATES (6)
         // Y1 = G
         for (int k = 0; k < NUMSTATES; k++) {
            // from i to k, take the highest probability and set to viterbi[t,k]
            max = 0;
            for (int i = 0; i < NUMSTATES; i++) {
               float prob = t_prob[i][k] * viterbi[t-1][i];
               if (prob > max)
                  max = prob;
            }
            viterbi[t][k] = max*e_prob[k][seqVal(t)];
         }
      }
      
      float max = 0;
      for (int i = 0; i < NUMSTATES; i++) {
         if (viterbi[size-1][i] > max) {
            max = viterbi[size-1][i];
            pathState[size-1] = i;
         }
      }
      
      probPath = max;
      
      for (int t = size-1; t > 0; t--) {
         // k = pathState[t]
         // t = t
         max = 0;
         for (int i = 0; i < NUMSTATES; i++) {
            float prob = viterbi[t-1][i]*t_prob[i][pathState[t]];
            if (prob > max) {
               max = prob;
               pathState[t-1] = i;
            }
         }
      }
   }
   
   // convert a nucleotide to it's integer value
   public static int seqVal(int loc) {
      char c = frag.charAt(loc);
      if (c == 'A')
         return A;
      else if (c == 'T')
         return T;
      else if (c == 'G')
         return G;
      else if (c == 'C')
         return C;
      System.err.println("Invalid nucleotide at location: " + (start + loc + 1));
      System.exit(-1);
      return -1;
   }
   
   // convert a state number to it's state name
   public static String stateVal(int val) {
      if (val == SS31)
         return "3'SS1";
      else if (val == SS32)
         return "3'SS2";
      else if (val == I)
         return "Intron";
      else if (val == SS51)
         return "5'SS1";
      else if (val == SS52)
         return "5'SS2";
      else if (val == E)
         return "Exon";
      System.err.println("Invalid state value of: " + val);
      System.exit(-1);
      return "";
   }
}
