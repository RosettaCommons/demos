/**
 * Created on Sep 1, 2009
 * Copyright Chris Poultney
 */
package parse;

import java.util.Arrays;
import java.util.Vector;
import util.Field;
import util.ResidueList;
import util.Stats;
import util.Field.FieldType;
import util.Residue;

/**
 * Produces median and variance as difference from wt.
 * @author crispy
 */
public class ParseV1 implements TsParser {
    boolean mean, median_compat, raw;
    FieldType[] inputMap = null;    // Handling instructions for input, indexed by field number
    FieldType[] outputMap = null;   // Instructions for output, which (may have) had fields added
    Vector<Double> wt = null;       // values from wt lines
    String pt_convert[] = { "nc", "ts", "lof", "?" };  // number-to-label conversion for v1-style -in.csv files
    boolean convertLabels;                 // true to convert numeric indices to labels, otherwise false
    public static String usage =
	"  V1: does median and variance calculation, output is median/var difference from wt\n" +
	"  -noconvertlabels: do not convert label indices to labels (v2 input: labels are already text)\n" +
	"  -mean: use mean instead of median\n" +
	"  -median_compat: calculate median the same way oocalc does for result compatibility\n" +
	"  -raw: output raw numbers instead of subtracting wt stats\n";
    /**
     * Construct a new V1 parser, which outputs median and difference in mutation/wt variance
     * for each run set.
     */
    public ParseV1() {
	mean = median_compat = raw = false;
	convertLabels = true;
    }
    /* (non-Javadoc)
     * @see parse.TsParser#registerOption(java.lang.String[], int)
     */
    @Override
    public int registerOption(String[] args, int idx) throws ArgumentException {
	int r;
	if(args[idx].equals("--usage") || args[idx].equals("--help")) {
	    throw new ArgumentException(usage);
	} else if(args[idx].equals("-mean")) {
	    mean = true;
	    r = idx + 1;
	} else if(args[idx].equals("-median_compat")) {
	    median_compat = true;
	    r = idx + 1;
	} else if(args[idx].equals("-raw")) {
	    raw = true;
	    r = idx + 1;
	} else if(args[idx].equals("-noconvertlabels")) {
	    convertLabels = false;
	    r = idx + 1;
	} else {
	    throw new ArgumentException("unknown v1 argument " + args[idx] + ", aborting\n" + usage);
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
		double f[] = new double[currmut.size()];
		for(int j = 0; j < currmut.size(); j++)
		    f[j] = Double.parseDouble(currmut.get(j)[i]);
		double s[] = Stats.getStats(f, median_compat);
		v.add(new Double(mean ? s[0] : s[1]));
		v.add(new Double(s[2]));
		break;
	    }
	}
	if(isWt)
	    wt = v;
	outputMutationLine(nat, mut, currmut, v, outputMap, wt);
    }
    /* (non-Javadoc)
     * @see parse.TsParser#finish()
     */
    @Override
    public void finish() {
	;  // nothing to do here
    }
    /**
     * Prints output line for a series of input lines pertaining to one mutation at one site.
     * @param nat native residue(s) at site
     * @param mut current residue(s) at site (could be native)
     * @param currmut each String in inner array is one field of one input line; one Vector entry per input line
     * @param processed array of values from {@link #processMutationLines(String, ResidueList, ResidueList, Vector, FieldType[])}
     * @param map field handling instructions
     * @param wt field values from wt line
     */
    public void outputMutationLine(ResidueList nat, ResidueList mut, Vector<String[]> currmut, Vector<Double> processed, FieldType[] map, Vector<Double> wt) {
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
		//			System.out.print(defFmt.format(processed[i+po]-wt[i+po]) + "," + defFmt.format(processed[i+po+1]-wt[i+po+1]));
		if(raw)
		    System.out.print(Field.defFmt.format(processed.get(i+po)) + "," + Field.defFmt.format(processed.get(i+po+1)));
		else
		    System.out.print(Field.defFmt.format(processed.get(i+po)-wt.get(i+po)) + "," + Field.defFmt.format(processed.get(i+po+1)-wt.get(i+po+1)));
		po++;
		break;
	    }
	}
	System.out.format("\n");
    }
    /**
     * Prints header line, including expansions for variance columns.
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
		System.out.format("%s,%s", header[i], header[i]+"V");
		break;
	    }
	}
	System.out.format("\n");
    }
}
