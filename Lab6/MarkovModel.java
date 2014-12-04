import java.util.*;
import java.io.*;

public class MarkovModel {
	static String fastaFile;
	static String gffFile;
	static int start;
	static int end;
	static float[][] transition;
	static float[][] emission;
	
	public static final int SS31 = 0;
	public static final int SS32 = 1;
	public static final int I = 2;
	public static final int SS51 = 3;
	public static final int SS52 = 4;
	public static final int E = 5;
	public static final int NumStates = 6;
	
	public static final int A = 0;
	public static final int T = 1;
	public static final int C = 2;
	public static final int G = 3;
	
   public static void main(String[] args) {
		transition = new float[6][6];
		emission = new float[6][4];
		parseConfigInput();
		int T = end - start;
		int k = NumStates;
		float[][] viterbi = new float[T][k];
		
		// TODO: lab 6
   }
   
  public static void parseConfigInput() {
	try {
		Scanner parser = new Scanner(new File("config.txt"));
		
		fastaFile = parser.nextLine().split("=")[1].trim();
		gffFile = parser.nextLine().split("=")[1].trim();
		start = Integer.parseInt(parser.nextLine().split("=")[1].trim());
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
}