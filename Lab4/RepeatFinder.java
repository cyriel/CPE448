import java.util.*;
import java.io.*;

public class RepeatFinder {
   static int start = 1, end = -1;

   public static ArrayList<String> crop(ArrayList<StrCntVst> list, int size) {
      ArrayList<String> ret = new ArrayList<String>();
      for (StrCntVst scv: list)
         if (scv.string.length() >= size && !ret.contains(scv.string.substring(0,size)))
            ret.add(scv.string.substring(0,size));
      return ret;
   }
   
   public static ArrayList<String> excludeN(ArrayList<String> list) {
      ArrayList<String> ret = new ArrayList<String>();
      for (String string: list)
         if (!string.contains("N"))
            ret.add(string);
      return ret;
   }
   
   public static double probability(String str) {
      double length = str.length();
      double prob = 1.0;
      int temp;
      temp = str.replace("A", "").length();
      if (temp != length)
         prob *= (length - temp)/length;
      temp = str.replace("T", "").length();
      if (temp != length)
         prob *= (length - temp)/length;
      temp = str.replace("G", "").length();
      if (temp != length)
         prob *= (length - temp)/length;
      temp = str.replace("C", "").length();
      if (temp != length)
         prob *= (length - temp)/length;
      return  prob;
   }
   
   public static double percentN(String str) {
      return 100.0*(str.length() - str.replace("A", "").length())/str.length();
   }
   
   public static void main2(String[] args) {
      SuffixTree2 tree = new SuffixTree2("GANANA", 3);
      ArrayList<StrCntVst> list = tree.getMaxRepStr();
      tree.printTree();
      ArrayList<String> shortened = crop(list, 3);
      ArrayList<StrCntVst> out = new ArrayList<StrCntVst>();
      for(String str: shortened)
         out.add(tree.search(str));
      
///*
      for(StrCntVst tuple: out) {
         System.out.print("String: " + tuple.string.replace(""+(char)0,"") + "\tCount: " + tuple.count + "\tLocations:");
         String loc = "";
         for (Integer i: tuple.visits)
            loc += ", " + i;
         System.out.println(loc.substring(1));
      }
//*/
   }
   
   public static void main(String[] args) {
      String config = "config.txt";
      Scanner parser;
      try {
         parser = new Scanner(new File(config));
         String fastafile = parser.nextLine();
         int sequenceSize = parser.nextInt();
         if (parser.hasNextInt()) {
            start = parser.nextInt();
            if (parser.hasNextInt())
               end = parser.nextInt();
            else
               start = 1;
         }
         start--;
         try {
            Scanner fasta = new Scanner(new File(fastafile));
            fasta.nextLine();
            String dna = "";
            while (fasta.hasNextLine())
               dna += fasta.nextLine();
            if (end == -1)
               end = dna.length();
            String frag = dna.substring(start,end);
            // this is where the magic needs to happen
            
            SuffixTree2 tree = new SuffixTree2(frag, sequenceSize);
//            tree.printTree();
            ArrayList<StrCntVst> list = tree.getMaxRepStr();
            for(StrCntVst tuple: list) {
//               System.out.print("String: " + tuple.string.replace(""+(char)0,"") + "\tCount: " + tuple.count + "\tLocations:");
               String loc = "";
               for (Integer i: tuple.visits)
                  loc += ", " + i;
//               System.out.println(loc.substring(1));
            }

            ArrayList<String> shortened = crop(list, sequenceSize);
            ArrayList<String> removedN = excludeN(shortened);
            ArrayList<StrCntVst> out = new ArrayList<StrCntVst>();
            for(String str: removedN) {
               out.add(tree.search(str));
            }
            
            Collections.sort(out);
            
            File txtOut = new File("output.txt");
            File csvOut = new File("output.csv");
            try {
               double percentN = percentN(frag);
               
               txtOut.createNewFile();
               csvOut.createNewFile();
               PrintWriter txtPrinter = new PrintWriter(new FileWriter(txtOut));
               PrintWriter csvPrinter = new PrintWriter(new FileWriter(csvOut));
               
               txtPrinter.println("Total %N = " + String.format("%.2f", percentN));
               csvPrinter.println("Total %N =," + percentN);
               txtPrinter.println("Total sequences checked = " + shortened.size());
               csvPrinter.println("Total sequences checked =," + shortened.size());
               txtPrinter.println("Total sequences checked without N = " + removedN.size());
               csvPrinter.println("Total sequences checked without N =," + removedN.size());
               
               txtPrinter.println(String.format("%20s%20s%20s","Sequence", "# of Repeats", "Probability"));
               csvPrinter.println("Sequence,# of Repeats,Probability");
               int ndx;
               for (ndx = 0; ndx < out.size(); ndx++) {
                  StrCntVst scv = out.get(ndx);
                  if (ndx < 10)
                     txtPrinter.println(String.format("%20s%20d%20.8f",scv.string, scv.count, probability(scv.string)));
                  csvPrinter.println(scv.string + "," + scv.count + "," + probability(scv.string));
               }
               
               txtPrinter.close();
               csvPrinter.close();
            } catch (IOException e) {
               System.err.println("IO error occured");
            }
         } catch (FileNotFoundException e) {
            System.err.println("FASTA file not found");
         }
      } catch (FileNotFoundException e) {
         System.err.println("Config file not found");
      }
      
   }
}

class SuffixTree2 {
   Node2 root;
/*
   public static void main(String[] args) {
      SuffixTree2 tree = new SuffixTree2("GANANA");
      ArrayList<StrCntVst> list = tree.getMaxRepStr();
      StrCntVst tuple = tree.search("ANAN");
      System.out.print("String: " + tuple.string + "\tCount: " + tuple.count + "\tLocations:");
      String loc = "";
      for (Integer i: tuple.visits)
         loc += ", " + i;
      System.out.println(loc.substring(1));
//      tree.printTree();

      for(StrCntVst tuple: list) {
         System.out.print("String: " + tuple.string + "\tCount: " + tuple.count + "\tLocations:");
         String loc = "";
         for (Integer i: tuple.visits)
            loc += ", " + i;
         System.out.println(loc.substring(1));
      }

   }
*/
   public SuffixTree2(String str, int size) {
      root = new Node2(true);
      for (int i = 0; i + size < str.length(); i++) {
//         if (i % 10000 == 0) System.out.println(i/10000);
         genST(root, str.substring(i, i+size), i);
      }
   }
   
   StrCntVst search(String str) {
      Node2 temp = root;
      String var = new String(str);
      int diff = temp.firstDifference(var);
      while (diff != var.length()) {
         var = var.substring(diff);
         if (temp.contains(var.charAt(0)))
            temp = temp.get(var.charAt(0));
         else
            return new StrCntVst(str, new ArrayList<Integer>());
         diff = temp.firstDifference(var);
      }
      return new StrCntVst(str, temp.visits);
   }
   
   ArrayList<StrCntVst> getMaxRepStr() {
      return root.getMaxRepStr();
   }
   
   void printTree() {
      printTree(0, root);
   }
   
   void printTabs(int n) {
      for (int i = 0; i < n; i++)
         System.out.print("\t");
   }
   
   void printTree(int n, Node2 node) {
      if (node.hasSuffixA) {
         printTabs(n);
         System.out.println(node.suffixA.root);
         printTree(n+1,node.suffixA);
      }
      if (node.hasSuffixT) {
         printTabs(n);
         System.out.println(node.suffixT.root);
         printTree(n+1,node.suffixT);
      }
      if (node.hasSuffixG) {
         printTabs(n);
         System.out.println(node.suffixG.root);
         printTree(n+1,node.suffixG);
      }
      if (node.hasSuffixC) {
         printTabs(n);
         System.out.println(node.suffixC.root);
         printTree(n+1,node.suffixC);
      }
      if (node.hasSuffixN) {
         printTabs(n);
         System.out.println(node.suffixN.root);
         printTree(n+1,node.suffixN);
      }
      if (node.hasEnd) {
         printTabs(n);
         System.out.println("X");
      }
   }
   
   void genST(Node2 node, String str, int start) {
      if (str.length() > 0) {
         // Node's visit will do: create nodes if necessary, shatter nodes if necessary, update the shattered child's visits with the node's visits, log visit to node and to child
         node.visit(str, start);
      }
   }
}

class Node2 {
//   Node2 prefix;
   String prefix = "";
   String root = "";
   Node2 suffixA, suffixT, suffixG, suffixC, suffixN;
   boolean hasChildren = false, hasSuffixA = false, hasSuffixT = false, hasSuffixG = false, hasSuffixC = false, hasSuffixN = false, hasEnd = true;
   boolean isRoot = false;
   int numChildren = 0;
   ArrayList<Integer> visits = new ArrayList<Integer>();
   
   public Node2(boolean isRoot) {
      this.isRoot = isRoot;
      hasEnd = false;
   }
   
   public Node2(Node2 node, String str) {
      prefix = node.prefix + node.root;
      root = str;
   }
   
   public Node2(Node2 node, String str, int visit) {
      prefix = node.prefix + node.root;
      root = str;
      visits.add(visit);
   }
   
   // GET MAX REPeating SubsTRing
   ArrayList<StrCntVst> getMaxRepStr() {
      ArrayList<StrCntVst> ret = new ArrayList<StrCntVst>();
      
//      if (!hasChildren)  // leaf node
//         return ret;

      // gets substring and count from each child
      if (hasSuffixA)
         ret.addAll(suffixA.getMaxRepStr());
      if (hasSuffixT)
         ret.addAll(suffixT.getMaxRepStr());
      if (hasSuffixG)
         ret.addAll(suffixG.getMaxRepStr());
      if (hasSuffixC)
         ret.addAll(suffixC.getMaxRepStr());
      if (hasSuffixN)
         ret.addAll(suffixN.getMaxRepStr());
      
      // if all children are leafs, they added nothing to the ArrayList
      // therefore the current node is not a leaf and is the longest substring here
      if (ret.isEmpty() && visits.size() > 1) { //&& numChildren > 1) {
         
         ret.add(new StrCntVst(prefix+root, visits));
      }
      
      return ret;
   }
   
/*   String getString() {
      Node2 temp = prefix;
      String ret = root;
      
      while (!temp.isRoot) {
         ret = temp.root + ret;
         temp = temp.prefix;
      }
      return ret;
   }*/
   
   int firstDifference(String str) {
      int i;
      for (i = 0; i < (root.length() > str.length() ? str.length() : root.length()); i++)
         if (root.charAt(i) != str.charAt(i))
            return i;
      return i;
   }
   
   // get loc from firstDifference()
   void shatter(int loc) {
      String temp = root.substring(loc);
      root = root.substring(0,loc);
      Node2 child = new Node2(this, temp);
      if (hasSuffixA) {
         child.addSuffixA(suffixA);
         hasSuffixA = false;
         suffixA = null;
      }
      if (hasSuffixT) {
         child.addSuffixT(suffixT);
         hasSuffixT = false;
         suffixT = null;
      }
      if (hasSuffixG) {
         child.addSuffixG(suffixG);
         hasSuffixG = false;
         suffixG = null;
      }
      if (hasSuffixC) {
         child.addSuffixC(suffixC);
         hasSuffixC = false;
         suffixC = null;
      }
      if (hasSuffixN) {
         child.addSuffixN(suffixN);
         hasSuffixN = false;
         suffixN = null;
      }
      hasEnd = false;
      child.visits = new ArrayList<Integer>(visits);
      child.numChildren = numChildren;
      
      char c = temp.charAt(0);
      if (c == 'A')
         addSuffixA(child);
      else if (c == 'T')
         addSuffixT(child);
      else if (c == 'G')
         addSuffixG(child);
      else if (c == 'C')
         addSuffixC(child);
      else if (c == 'N')
         addSuffixN(child);
      numChildren = 1;
   }
   
   // unsafe adding of suffixes. Tread carefully
   private void addSuffixA(Node2 node) {
      suffixA = node;
      hasChildren = hasSuffixA = true;
   }
   
   private void addSuffixT(Node2 node) {
      suffixT = node;
      hasChildren = hasSuffixT = true;
   }
   
   private void addSuffixG(Node2 node) {
      suffixG = node;
      hasChildren = hasSuffixG = true;
   }
   
   private void addSuffixC(Node2 node) {
      suffixC = node;
      hasChildren = hasSuffixC = true;
   }
   
   private void addSuffixN(Node2 node) {
      suffixN = node;
      hasChildren = hasSuffixN = true;
   }
   
   void visit(String str, int start) {
         int diff = firstDifference(str);
         if (root.length() > diff) {
            shatter(diff);
         }
         if (!isRoot) {
            visits.add(start);
         }
            
         if (str.length() == diff) {
            hasEnd = true;
            numChildren++;
            return;
         }
         str = str.substring(diff);
      char c = str.charAt(0);
      
      if (c == 'A') {
         if (!hasSuffixA) {
            suffixA = new Node2(this, str, start);
            hasChildren = hasSuffixA = true;
            numChildren++;
         } else // painful creation here
            suffixA.visit(str, start);
      } else if (c == 'T') {
         if (!hasSuffixT) {
            suffixT = new Node2(this, str, start);
            hasChildren = hasSuffixT = true;
            numChildren++;
         } else // painful creation here
            suffixT.visit(str, start);
      } else if (c == 'G') {
         if (!hasSuffixG) {
            suffixG = new Node2(this, str, start);
            hasChildren = hasSuffixG = true;
            numChildren++;
         } else // painful creation here
            suffixG.visit(str, start);
      } else if (c == 'C') {
         if (!hasSuffixC) {
            suffixC = new Node2(this, str, start);
            hasChildren = hasSuffixC = true;
            numChildren++;
         } else // painful creation here
            suffixC.visit(str, start);
      } else if (c == 'N') {
         if (!hasSuffixN) {
            suffixN = new Node2(this, str, start);
            hasChildren = hasSuffixN = true;
            numChildren++;
         } else // painful creation here
            suffixN.visit(str, start);
      } else
         System.err.println("Error: Invalid character input.");
   }
   
   boolean contains(char c) {
      if (c == 'A')
         return hasSuffixA;
      else if (c == 'T')
         return hasSuffixT;
      else if (c == 'G')
         return hasSuffixG;
      else if (c == 'C')
         return hasSuffixC;
      else if (c == 'N')
         return hasSuffixN;
      else
         System.err.println("Error: Invalid character input.");
      return false;
   }
   
   Node2 get(char c) {
      if (c == 'A')
         return suffixA;
      else if (c == 'T')
         return suffixT;
      else if (c == 'G')
         return suffixG;
      else if (c == 'C')
         return suffixC;
      else if (c == 'N')
         return suffixN;
      else
         System.err.println("Error: Invalid character input.");
      return null;
   }
}

final class StrCntVst implements Comparable<StrCntVst>{
   final String string;
   final int count;
   final ArrayList<Integer> visits;
   
   public StrCntVst(String string, ArrayList<Integer> visits) {
      this.string = string;
      this.count = visits.size();
      this.visits = visits;
   }
   
   public int compareTo(StrCntVst other) {
      return other.count - count;
   }
}
