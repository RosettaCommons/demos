/**
 * Created on Mar 7, 2010
 * Copyright Chris Poultney
 */
package app;



import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.net.URL;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Pattern;
import util.PSSM;

/**
 * @author crispy
 *
 */
public class Test {
    static int INFO = 1, SCORE = 2, PERCENT = 4;
    public static void main(String[] args) {
//	testPatterns();
//	testMap();
	try {
	    testPSSM(SCORE);
//	    doThing();
	} catch(Exception e) {
	    e.printStackTrace();
	    System.exit(1);
	}
    }
    public static void testMap() {
	TreeMap<Integer,Integer> corresp = new TreeMap<Integer,Integer>();
	corresp.put(new Integer(38), new Integer(1));
	corresp.put(new Integer(224), new Integer(98));
	int[] from = { 38, 134, 224, 302 };
	for(int i = 0; i < from.length; i++) {
	    Map.Entry<Integer,Integer> m = corresp.floorEntry(new Integer(from[i]));
	    int offset = m.getKey().intValue() - m.getValue().intValue();
	    System.out.println(from[i] + " ==> " + (from[i] - offset));
	}
    }
    public static void testPatterns() {
	String test = "        { 433430569842969, 10, -17 },";
	String[] f = test.split("[\\s{},]+");
	for(int i = 0; i < f.length; i++)
	    System.out.format("%d: %s\n", i, f[i]);
	Pattern stripP = Pattern.compile("^\\s+(.*[^,]),*$");
	String test2 = "    isProtein TRUE,";
	System.out.println(stripP.matcher(test2).replaceAll("$1"));
    }
    public static void testPSSM(int out) throws IOException {
	PSSM pssm = PSSM.parsePSSM(new URL("file:../../../yeast2/pssm/YFL039C.cp.gz"));
	System.out.println("success");
	// show info content
	if((out & INFO) > 0) {
	    for(int i = 0; i < pssm.getCols(); i++)
		System.out.println(i + ": " + pssm.getInfoContent(i));
	}
	// pre-processing for matrix display
	String order = "ARNDCQEGHILKMFPSTWYV";
	int[] corr = new int[20];
	for(int i = 0; i < corr.length; i++)
	    corr[i] = PSSM.STDAA.indexOf(order.charAt(i));
	// show score matrix
	if((out & SCORE) > 0) {
	    System.out.print("\t");
	    for(int i = 0; i < order.length(); i++)
		System.out.format("%4s", order.charAt(i));
	    System.out.println();
	    for(int c = 1; c <= pssm.getCols(); c++) {
		System.out.format("%4d %c\t", c, pssm.getSequenceAt(c));
		for(int i = 0; i < order.length(); i++)
		    System.out.format("%4d", pssm.getScore(c, order.charAt(i)));
		System.out.println();
	    }
	}
	// show percentages
	if((out & PERCENT) > 0) {
	    for(int i = 0; i < order.length(); i++)
		System.out.format("%4s", order.charAt(i));
	    System.out.println();
	    for(int c = 0; c < pssm.getCols(); c++) {
		for(int i = 0; i < order.length(); i++) {
		    System.out.format("%4d", (int)Math.round(pssm.getWeightedFreq(c, order.charAt(i)) * 100));
		}
		System.out.println();
	    }
	}
    }
    public static void doThing() throws IOException {
	BufferedReader in = null;
	try {
	    in = new BufferedReader(new FileReader("tmp.txt"));
	    String s;
	    while((s = in.readLine()) != null)
		System.out.println(s);
	} catch(IOException ie) {
	    throw ie;
	} finally {
	    if(in != null)
		in.close();
	    System.err.println("closed");
	}
    }
}
