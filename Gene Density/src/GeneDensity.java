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
		return rec.direction ? dna.substring(rec.start,rec.stop) : reverseComplement(dna.substring(rec.stop,rec.start));
	}
	
	public static String extractGeneDNA(String dna, Isoform iso) {
		return iso.direction ? dna.substring(iso.start, iso.stop) : reverseComplement(dna.substring(iso.stop, iso.start));
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
	public int isNested(Isoform iso1, Isoform iso2) {
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
	public int[] nested(Isoform iso1, Isoform iso2, int hint) {
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
	      range = (iso1.direction) ? new int[] {iso1.start,iso1.stop} : new int[] {iso1.stop, iso1.start};
      else if (hint == 7)
	      range = new int[] {iso2.start,iso1.stop};
      else if (hint == 8)
	      range = new int[] {iso1.start,iso2.start};
	   else if (hint == 10)
	      range = new int[] {iso2.stop,iso1.stop};
	   else if (hint == 11)
	      range = new int[] {iso1.start,iso2.stop};
	   else if (hint == 9 || hint == 12)
	      range = (iso2.direction) ? new int[] {iso2.start,iso2.stop} : new int[] {iso2.stop, iso2.start};
	   else
	      range = new int[] {0,0};
	      
      return range;
	}
	
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
