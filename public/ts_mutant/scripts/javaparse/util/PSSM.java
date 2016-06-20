/**
 * Created on Mar 7, 2010
 * Copyright Chris Poultney
 */
package util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.net.URL;
import java.util.Hashtable;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 * @author crispy
 *
 */
public class PSSM {
    public static final String STDAA = "-ABCDEFGHIKLMNPQRSTVWXYZU*JO";
    public static final int ROWS = STDAA.length();
    int cols;
    double[][] weightedFreqs;
    double[] infoContent;
    int[][] scores;
    String sequence;
    Hashtable<String,Integer> aaIndex;
    private PSSM(int cols, double[][] wf, double[] ic, int[][] s, String seq) {
	this.cols = cols;
	this.weightedFreqs = wf;
	this.infoContent = ic;
	this.scores = s;
	this.sequence = seq;
	aaIndex = new Hashtable<String,Integer>();
	for(int i = 0; i < STDAA.length(); i++)
	    aaIndex.put(STDAA.charAt(i)+"", new Integer(i));
    }
    public int getCols() { return cols; }
    public int getScore(int c, String aa) { return scores[c-1][aaIndex.get(aa).intValue()]; };
    public int getScore(int c, char aa) { return scores[c-1][aaIndex.get(aa+"").intValue()]; };
    public int getScore(int c, int r) { return scores[c-1][r]; }
    public double getWeightedFreq(int c, String aa) { return weightedFreqs[c-1][aaIndex.get(aa).intValue()]; }
    public double getWeightedFreq(int c, char aa) { return weightedFreqs[c-1][aaIndex.get(aa+"").intValue()]; }
    public double getWeightedFreq(int c, int r) { return weightedFreqs[c-1][r]; }
    public double getInfoContent(int c) { return infoContent[c-1]; }
    public String getSequence() { return sequence; }
    public char getSequenceAt(int c) { return sequence.charAt(c-1); }
    public static PSSM parsePSSM(URL src) throws IOException {
	PSSMReader data = new PSSMReader(src);
	return new PSSM(data.cols, data.weightedFreqs, data.infoContent, data.scores, data.sequence);
    }
    static class PSSMReader {
	URL src;
	int cols;
	boolean byRow;
	double[][] weightedFreqs;
	double[] infoContent;
	int[][] scores;
	String sequence;
	static final Pattern stripP = Pattern.compile("^\\s+(.*[^,]),*$");
	public PSSMReader() {
	    src = null;
	    cols = 0;
	    weightedFreqs = null;
	    infoContent = null;
	    scores = null;
	    sequence = null;
	}
	public PSSMReader(URL src) throws IOException {
	    this.src = src;
	    LineNumberReader in = new LineNumberReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(src.openStream()))));
	    try {
		parse(in);
	    } catch(IOException e) {
		throw new IOException("Error at line " + in.getLineNumber() +": ", e);
	    } finally {
		in.close();
	    }
	}
	public void parse(LineNumberReader in) throws IOException {
	    String s = in.readLine();
	    if(!s.startsWith("PssmWithParameters"))
		throw new IOException("URL does not appear to describe ascii_pssm file with parameters: " + src);
	    while((s = in.readLine()) != null) {
		String[] f = parseLine(s);
		if(f[0].equals("pssm")) {
		    parsePssm(in);
		} else if(f[0].equals("params")) {
		    readUntilBlockEnd(in);
		} else if(f[0].equals("}")) {
		    break; // end of section
		} else {
		    throw new IOException("Unexpected main block parameter: \"" + f[0] + "\"");
		}
	    }
	    in.close();
	    if(weightedFreqs== null)
		throw new IOException("Missing parameter: weightedFreqs");
	    if(infoContent == null)
		throw new IOException("Missing parameter: infoContent");
	    if(scores == null)
		throw new IOException("Missing parameter: scores");
	    if(sequence == null)
		throw new IOException("Missing parameter: sequence");
	}
	public void parsePssm(LineNumberReader in) throws IOException {
	    String s;
	    String last = null;
	    while((s = in.readLine()) != null) {
		String[] f = parseLine(s);
		if(f[0].equals("isProtein")) {
		    boolean b = Boolean.parseBoolean(f[1]);
		    if(!b)
			throw new IOException("isProtein flag is false: " + f[1]);
		} else if(f[0].equals("numRows")) {
		    int i = Integer.parseInt(f[1]);
		    if(i != ROWS)
			throw new IOException("Unexpected number of rows: " + i);
		} else if(f[0].equals("numColumns")) {
		    cols = Integer.parseInt(f[1]);
		} else if(f[0].equals("byRow")) {
		    boolean b = Boolean.parseBoolean(f[1]);
		    if(b)
			throw new IOException("Cannot handle files sorted by row");
		} else if(f[0].equals("query")) {
		    parseQuery(in);
		} else if(f[0].equals("intermediateData")) {
		    parseIntermediateData(in);
		} else if(f[0].equals("finalData")) {
		    parseFinalData(in);
		} else if(f[0].equals("identifier")) {
		    ;
		} else if(f[0].equals("}")) {
		    break; // end of section
		} else {
		    throw new IOException("Unexpected pssm block parameter: " + f[0] + ", last: " + last);
		}
		last = s;
	    }
	}
	public void parseQuery(LineNumberReader in) throws IOException {
	    String s;
	    while((s = in.readLine()) != null) {
		String[] f = parseLine(s);
		if(f[0].equals("id")) {
		    readUntilBlockEnd(in);
		} else if(f[0].equals("inst")) {
		    parseInst(in);
		} else if(f[0].equals("descr")) {
		    readUntilBlockEnd(in);
		} else if(f[0].equals("}")) {
		    break; // end of section
		} else {
		    throw new IOException("Unexpected query block parameter: " + f[0]);
		}
	    }
	}
	public void parseInst(LineNumberReader in) throws IOException {
	    String s;
	    while((s = in.readLine()) != null) {
		String[] f = parseLine(s);
		if(f[0].equals("repr")) {
		    continue;
		} else if(f[0].equals("mol")) {
		    continue;
		} else if(f[0].equals("length")) {
		    int l = Integer.parseInt(f[1]);
		    if(l != cols)
			throw new IOException("advertised sequence length does not equal number of columns: " + l + " != " + cols);
		} else if(f[0].equals("seq-data")) {
		    // this is different because it's multi-line data delimited with single quotes instead of brackets
		    // first line format is seq-data ncbistdaa '....
		    // query lines continue until pattern ...'H
		    String seq = f[2];
		    while(!seq.endsWith("'H"))
			seq += in.readLine();
		    // strip single quotes and terminating 'H'
		    Pattern pat = Pattern.compile("'(.*)'H");
		    seq = pat.matcher(seq).replaceAll("$1");
		    StringBuffer stmp = new StringBuffer();
		    for(int i = 0; i < seq.length(); i+=2) {
			int v = Integer.parseInt(seq.substring(i, i+2), 16);
			stmp.append(STDAA.charAt(v));
		    }
		    sequence = stmp.toString();
		    if(sequence.length() != cols)
			throw new IOException("actual sequence length does not equal number of columns: " + sequence.length() + " != " + cols);
		} else if(f[0].equals("}")) {
		    break;
		} else {
		    throw new IOException("Unexpected inst block paramater: " + f[0]);
		}
	    }
	}
	public void parseIntermediateData(LineNumberReader in) throws IOException {
	    String s;
	    while((s = in.readLine()) != null) {
		String[] f = parseLine(s);
		if(f[0].equals("weightedResFreqsPerPos")) {
		    parseWeightedFreqs(in);
		} else if(f[0].equals("informationContent")) {
		    parseInfoContent(in);
		} else if(f[0].equals("freqRatios")) {
		    readUntilBlockEnd(in);
		} else if(f[0].equals("gaplessColumnWeights")) {
		    readUntilBlockEnd(in);
		} else if(f[0].equals("sigma")) {
		    readUntilBlockEnd(in);
		} else if(f[0].equals("intervalSizes")) {
		    readUntilBlockEnd(in);
		} else if(f[0].equals("numMatchingSeqs")) {
		    readUntilBlockEnd(in);
		} else if(f[0].equals("}")) {
		    break; // end of section
		} else {
		    throw new IOException("Unexpected intermediateData block parameter: " + f[0]);
		}
	    }
	}
	public void parseFinalData(LineNumberReader in) throws IOException {
	    String s;
	    while((s = in.readLine()) != null) {
		String[] f = parseLine(s);
		if(f[0].equals("scores")) {
		    parseScores(in);
		} else if(f[0].equals("lambda")) {
		    ;  // single line name/value(REAL)
		} else if(f[0].equals("kappa")) {
		    ;  // single line name/value(REAL)
		} else if(f[0].equals("h")) {
		    ;  // single line name/value(REAL)
		} else if(f[0].equals("lambdaUngapped")) {
		    ;  // single line name/value(REAL)
		} else if(f[0].equals("kappaUngapped")) {
		    ;  // single line name/value(REAL)
		} else if(f[0].equals("hUngapped")) {
		    ;  // single line name/value(REAL)
		} else if(f[0].equals("}")) {
		    break; // end of section
		} else {
		    throw new IOException("Unexpected finalData block parameter: " + f[0]);
		}
	    }
	}
	public void parseWeightedFreqs(LineNumberReader in) throws IOException {
	    weightedFreqs = new double[cols][ROWS];
	    String s;
	    int cnt = 0;
	    while((s = in.readLine()) != null) {
		String[] f = parseLine(s);
		if(f[0].equals("}")) {
		    break;
		} else { // assume usual triple
		    try {
			weightedFreqs[cnt/ROWS][cnt%ROWS] = parseReal(s);
			cnt++;
		    } catch(Exception e) {
			throw new IOException(e);
		    }
		}
	    }
	    if(cnt != cols * ROWS)
		throw new IOException("Unexpected number of weighted frequency entries: got " + cnt + ", expected " + (cols*ROWS));
	}
	public void parseInfoContent(LineNumberReader in) throws IOException {
	    infoContent = new double[cols];
	    String s;
	    int cnt = 0;
	    while((s = in.readLine()) != null) {
		String[] f = parseLine(s);
		if(f[0].equals("}")) {
		    break;
		} else { // assume usual triple
		    try {
			infoContent[cnt] = parseReal(s);
			cnt++;
		    } catch(Exception e) {
			throw new IOException(e);
		    }
		}
	    }
	    if(cnt != cols)
		throw new IOException("Unexpected number of information content entries: got " + cnt + ", exptected " + cols);
	}
	public void parseScores(LineNumberReader in) throws IOException {
	    scores = new int[cols][ROWS];
	    String s;
	    int cnt = 0;
	    while((s = in.readLine()) != null) {
		String[] f = parseLine(s);
		if(f[0].equals("}")) {
		    break;
		} else { // assume INTEGER
		    try {
			scores[cnt/ROWS][cnt%ROWS] = Integer.parseInt(f[0]);
			cnt++;
		    } catch(Exception e) {
			throw new IOException(e);
		    }
		}
	    }
	    if(cnt != cols * ROWS)
		throw new IOException("Unexpected number of score entries: got " + cnt + ", expected " + (cols*ROWS));
	}
	public void readUntilBlockEnd(LineNumberReader in) throws IOException {
	    int indent = 1;
	    String s;
	    while(indent > 0 && (s = in.readLine()) != null) {
		String[] f = parseLine(s);
		for(int i = 0; i < f.length; i++)
		    if(f[i].equals("{"))
			indent++;
		    else if(f[i].equals("}"))
			indent--;
	    }
	}
	public String[] parseLine(String s) {
	    // types:
	    //   name {
	    //   name value[,]
	    //   { tuple }[,]
	    //   }
//	    System.err.println("before: " + s);
	    String p = stripP.matcher(s).replaceAll("$1");
//	    System.err.println(" after: " + p);
	    if(p.startsWith("{"))
		return new String[] { p };
	    else
		return p.split("\\s");
	}
	public double parseReal(String s) {
	    // format is "{ mantissa, base, exp }"
	    //  f[0] is blank because of leading separators, f[1]-f[3] contain data
	    String[] f = s.split("[\\s{},]+");
	    return Double.parseDouble(f[1]) * Math.pow(Double.parseDouble(f[2]), Double.parseDouble(f[3]));
	}
    }
}
