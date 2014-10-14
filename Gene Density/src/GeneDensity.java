import java.util.ArrayList;

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
	public static String extractRecordDNA(String dna, Record rec) {
		String extracted = null;
		//some logic happens here
		return rec.direction ? extracted : reverseComplement(extracted);
	}
	
	public static String extractGeneDNA(String dna, Isoform iso) {
		String extracted = null;
		if 
		//some logic happens here
		return iso.direction ? extracted : reverseComplement(extracted);
	}
	
	public static String reverseComplement(String dna) {
		char[] reverse = new StringBuilder(dna).reverse().toString().toCharArray();
		for (int i = 0; i<reverse.length; i++) {
			if (reverse[i] == 'A')
				reverse[i] = 'T';
			else if (reverse[i] == 'T')
				reverse[i] = 'A';
			else if (reverse[i] == 'G')
				reverse[i] = 'C';
			else if (reverse[i] == 'C')
				reverse[i] = 'G';
		}
		return reverse.toString();
	}
	
//	nested genes: overlap start and stop
//	intergenic region: sort list by smallest index (either start or stop)
}

class Isoform {
	ArrayList<Record> records = new ArrayList<Record>();
	int start;
	int stop;
	boolean direction;
	
	public void setStart(int start) {
		this.start = start;
	}
	
	public void setStop(int stop) {
		this.stop = stop;
	}
	
	public void setDirection(boolean Direction) {
		this.direction = direction;
	}
}

class Record {
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