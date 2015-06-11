/**
 * Created on Sep 1, 2009
 * Copyright Chris Poultney
 */
package parse;

import java.util.Arrays;
import java.util.Vector;
import util.Field;
import util.Residue;
import util.ResidueList;
import util.Stats;
import util.Field.FieldType;

/**
 * Produces 3 output terms per input terms: position of 1st-3rd quartiles of a mutation
 * within the wt distribution, normalized to 0-1 range.
 * @author crispy
 * @deprecated replaced by {@link ParseV2} (formerly ParseV3), which produces the same output
 * when run with the -local flag as this class (formerly ParseV2) did.
 */
public class ParseV2Old implements TsParser {
    FieldType[] inputMap = null;    // Handling instructions for input, indexed by field number
    FieldType[] outputMap = null;   // Instructions for output, which (may have) had fields added
    Vector<double[]> wt;            // values from wt lines
    String pt_convert[] = { "nc", "ts", "lof", "?" };  // number-to-label conversion for v1-style -in.csv files
    boolean convertLabels;          // true to convert numeric indices to labels, otherwise false
    public static String usage =
	"  V2: does quartile calculation, outputs position of 1st-3rd quartiles of mutation relative to wt\n" +
	"  -convertlabels: convert label indices to labels (v1 input: labels are indices, not text)\n";
    /**
     * Construct a new V2 parser, which outputs median and difference in mutation/wt variance
     * for each run set.
     * @param args argument string passed to executable
     * @param idx index of first argument applicable to this class
     * @throws ArgumentException
     */
    public ParseV2Old() {
	wt = new Vector<double[]>();
	convertLabels = false;
    }
    /* (non-Javadoc)
     * @see parse.TsParser#registerOption(java.lang.String[], int)
     */
    @Override
    public int registerOption(String[] args, int idx) throws ArgumentException {
	int r;
	if(args[idx].equals("--usage") || args[idx].equals("--help")) {
	    throw new ArgumentException(usage);
	} else if(args[idx].equals("-noconvertlabels")) {
	    convertLabels = false;
	    r = idx + 1;
	} else {
	    throw new ArgumentException("unknown v2 argument " + args[idx] + ", aborting\n" + usage);
	}
	return r;
    }
    @Override
    public String getExtraFields() {
	return "aminochange,aminochangeglenn";
    }
    /* (non-Javadoc)
     * @see parse.TsParser#processHeader(java.lang.String[], util.Field.FieldType[])
     */
    @Override
    public void processHeader(String[] header, FieldType[] inputMap) {
	this.inputMap = inputMap;
	outputMap = Field.getFieldMap(header);
	outputHeaderLine(header, outputMap);
	for(int i = 0; i < outputMap.length; i++)
	    wt.add(new double[0]);
    }

    /* (non-Javadoc)
     * @see parse.TsParser#processLine(java.lang.String, util.ResidueList, java.lang.String[])
     */
    @Override
    public void processLine(String name, ResidueList mut, String[] fields) {
	if(convertLabels)
	    fields[fields.length-1] = pt_convert[Integer.parseInt(fields[fields.length-1])];
    }
    /* (non-Javadoc)
     * @see parse.TsParser#processMutation(java.lang.String, util.ResidueList, util.ResidueList, java.util.Vector, boolean)
     */
    @Override
    public void processMutation(String name, ResidueList mut, ResidueList nat, Vector<String[]> currmut, boolean isWt, String cwd, boolean useCorresp) {
	Vector<Double> v = new Vector<Double>();
	for(int i = 0; i < inputMap.length; i++) {
	    switch(inputMap[i]) {
	    case IGNORE:
	    case DONTCONVERT:
		v.add(new Double(0));
		break;
	    case CONVERT:
		// parse strings into array
		double f[] = new double[currmut.size()];
		for(int j = 0; j < currmut.size(); j++)
		    f[j] = Double.parseDouble(currmut.get(j)[i]);
		Arrays.sort(f);
		if(isWt)
		    wt.set(i, f);
		double q1 = Stats.pctlerp(.25, f), q2 = Stats.pctlerp(.50, f), q3 = Stats.pctlerp(.75, f);
		v.add(new Double(Stats.normfind(q1, wt.get(i))));
		v.add(new Double(Stats.normfind(q2, wt.get(i))));
		v.add(new Double(Stats.normfind(q3, wt.get(i))));
		break;
	    }
	}
	outputMutationLine(nat, mut, currmut, v, outputMap);
    }
    /* (non-Javadoc)
     * @see parse.TsParser#finish()
     */
    @Override
    public void finish() {
	;	// nothing to do here
    }
    /**
     * Prints header line, including expansions for quartile columns.
     * @param header separated fields of header line
     * @param map field handling instructions
     */
    public void outputHeaderLine(String[] header, FieldType[] map) {
	for(int i = 0; i < map.length; i++) {
	    switch(map[i]) {
	    case IGNORE: break;
	    case DONTCONVERT:
		if(i > 0) System.out.format(",");
		System.out.format("%s", header[i]);
		break;
	    case CONVERT:
		if(i > 0) System.out.format(",");
		System.out.format("%s,%s,%s", header[i]+"Q1", header[i]+"Q2", header[i]+"Q3");
		break;
	    }
	}
	System.out.format("\n");
    }
    /**
     * Prints output line for a series of input lines pertaining to one mutation at one site.
     * @param nat native residue(s) at site
     * @param mut current residue(s) at site (could be native)
     * @param currmut each String in inner array is one field of one input line; one Vector entry per input line
     * @param processed array of values from {@link #processMutationLines(String, ResidueList, ResidueList, Vector, FieldType[])}
     * @param map field handling instructions
     */
    public void outputMutationLine(ResidueList nat, ResidueList mut, Vector<String[]> currmut, Vector<Double> processed, FieldType[] map) {
	// this creates the representative line for DONTCONVERT fields
	// any per-mutation fields added in header processing should be inserted here
	// must be coordinated with field addition in header processing
	Vector<String> def = new Vector<String>(Arrays.asList(currmut.get(0)));
	double ac = Residue.getAminoChange(nat, mut);
	double acg = Residue.getAminoChangeGlenn(nat, mut);
	def.add(3, Field.chgFmt.format(ac));
	def.add(4, Field.chgFmt.format(acg));
	processed.add(3, new Double(0));
	processed.add(4, new Double(0));
	int po = 0;  // offset between indices for processed and currmut
	for(int i = 0; i < map.length; i++) {
	    switch(map[i]) {
	    case IGNORE: break;
	    case DONTCONVERT:
		if(i > 0) System.out.format(",");
		System.out.format("%s", def.get(i));
		break;
	    case CONVERT:
		if(i > 0) System.out.format(",");
		System.out.print(Field.defFmt.format(processed.get(i+po))
			+ "," + Field.defFmt.format(processed.get(i+po+1))
			+ "," + Field.defFmt.format(processed.get(i+po+2)));
		po+=2;
		break;
	    }
	}
	System.out.format("\n");
    }
}
