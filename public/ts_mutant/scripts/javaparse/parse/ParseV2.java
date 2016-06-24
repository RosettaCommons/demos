/**
 * Created on Sep 16, 2009
 * Copyright Chris Poultney
 */
package parse;


import java.util.Arrays;
import java.util.TreeSet;
import java.util.Vector;
import util.Field;
import util.Residue;
import util.ResidueList;
import util.Sequence;
import util.Stats;
import util.Field.FieldType;
import util.Sequence.SeqField;
import util.Sequence.SequenceAnalyzer;
import util.Sequence.SequenceException;

/**
 * Localrelax output parsing.
 * <p>Behavior of this class is by design a little odd.
 * Processing is done in two parts to accommodate the global option.  In both
 * local and global modes, input is split into String ({@link FieldType#DONTCONVERT}) and
 * numeric ({@link FieldType#CONVERT}) vectors by
 * {@link #processMutation(String, ResidueList, ResidueList, Vector, boolean)},
 * and values for the quartiles are stored as {@link GlobalEntry} objects.  Calculation
 * of the percentile values w.r.t. wt of the mutation quartiles is completed in
 * {@link #outputSemiprocessedMutationLine(GlobalEntry, Vector, FieldType[])}.</p>
 * <p>In local mode, line output is done directly from
 * {@link #processMutation(String, ResidueList, ResidueList, Vector, boolean)},
 * using the most recent {@link GlobalEntry} object and wt lines stored in
 * {@link #wt_local}.  In global mode, all
 * wt lines and the {@link GlobalEntry} objects that hold partially completed calculations
 * are stored in {@link #wt_global} and {@link #semiProcessed}, respectively, and final
 * calculation and output are done from the {@link #finish()} method after input has
 * been exhausted.</p>
 * 3/1/2010: Global mode has been disabled until it can be fixed.
 * Global mode works as designed, but is semantically flawed. Collecting wt lines
 * over all positions works as designed, but is not correct for more than one protein.
 * Global mode should be eliminated or fixed to collect per-protein instead of over the
 * entire input. What was intended to be global processing works in local mode because the
 * wt lines are repeated for each input position even though they are the same across all
 * positions of a protein for global runs.
 * @author crispy
 */
public class ParseV2 implements TsParser {
    FieldType[] inputMap = null;    // Handling instructions for input, indexed by field number
    FieldType[] outputMap = null;   // Instructions for output, which (may have) had fields added
    enum RunType { LOCAL, GLOBAL };
    RunType type;	            // Type of V2 processing being done
    Vector<double[]> wt_local;      // values from wt lines
    Vector<String[]> wt_global;     // raw fields for all wt lines
    Vector<GlobalEntry> semiProcessed;
    SeqField[] seqopt = { SeqField.B, SeqField.P, SeqField.Q, SeqField.R, SeqField.S }; // options for sequence-based term processing and output
    SequenceAnalyzer seqan;
    static String SUBDIR1 = "pssm2";
    static String SUBDIR2 = "pssm";
    TreeSet<String> badProt;
    public static String usage =
	"  V2: does quartile calculation, outputs position of 1st-3rd quartiles of mutation relative to wt\n" +
	"  -local: (default) uses only relax runs of the native at that position for wt statistics\n" +
	"  -global: uses relax runs of native structure across all positions for wt statistics\n" +
    	"  -s: sub-option for sequence-based terms (default: BPQRS)\n" +
    	"    terms as follows:\n";
    static {
	SeqField[] f = SeqField.values();
	for(int i = 0; i < f.length; i++)
	    usage += "    " + f[i].name() + ": " + f[i].getFieldName() + "\n";
    }
    public ParseV2() {
	type = RunType.LOCAL;
	wt_local = new Vector<double[]>();
	wt_global = new Vector<String[]>();
	semiProcessed = new Vector<GlobalEntry>();
	seqan = new SequenceAnalyzer();
	badProt = new TreeSet<String>();
    }
    /* (non-Javadoc)
     * @see parse.TsParser#registerOption(java.lang.String[], int)
     */
    @Override
    public int registerOption(String[] args, int idx) throws ArgumentException {
	int r;
	if(args[idx].equals("--usage") || args[idx].equals("--help")) {
	    throw new ArgumentException(usage);
	} else if(args[idx].equals("-local")) {
	    type = RunType.LOCAL;
	    r = idx + 1;
	} else if(args[idx].equals("-global")) {
	    throw new ArgumentException("-global flag has been disabled, aborting");
//	    type = RunType.GLOBAL;
//	    r = idx + 1;
	} else if(args[idx].equals("-s")) {
	    Vector<SeqField> v = new Vector<SeqField>();
	    String s = args[idx+1];
	    for(int i = 0; i < s.length(); i++) {
		try {
		    v.add(Sequence.getField(s.charAt(i)));
		} catch(IllegalArgumentException ia) {
		    throw new ArgumentException("invalid -s option \"" + s.charAt(i) + "\", aborting\n" + usage);
		}
	    }
	    seqopt = (SeqField[])v.toArray(new SeqField[0]);
	    r = idx + 2;
	} else {
	    throw new ArgumentException("unknown v1 argument " + args[idx] + ", aborting\n" + usage);
	}
	return r;
    }
    @Override
    public String getExtraFields() {
	String s = "";
	for(int i = 0; i < seqopt.length; i++)
	    s += (s == "" ? "" : ",") + seqopt[i].getFieldName();
	return s;
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
	    wt_local.add(new double[0]);
    }
    /* (non-Javadoc)
     * @see parse.TsParser#processLine(java.lang.String, util.ResidueList, java.lang.String[])
     */
    @Override
    public void processLine(String name, ResidueList mut, String[] fields) {
	;	// this method intentionally left blank
    }
    /* (non-Javadoc)
     * @see parse.TsParser#processMutation(java.lang.String, util.ResidueList, util.ResidueList, java.util.Vector, boolean)
     */
    @Override
    public void processMutation(String name, ResidueList mut, ResidueList nat, Vector<String[]> currmut, boolean isWt, String cwd, boolean useCorresp) {
	if(isWt)
	    if(type == RunType.LOCAL) {
		wt_local = processWT(currmut);
	    } else {
		wt_global.addAll(currmut);
	    }
	Vector<String> vs = new Vector<String>();
	Vector<Double> vd = new Vector<Double>();
	for(int i = 0; i < inputMap.length; i++) {
	    switch(inputMap[i]) {
	    case IGNORE: break;
	    case DONTCONVERT: vs.add(currmut.get(0)[i]); break;
	    case CONVERT:
		double f[] = new double[currmut.size()];
		for(int j = 0; j < currmut.size(); j++)
		    f[j] = Double.parseDouble(currmut.get(j)[i]);
		Arrays.sort(f);
		vd.add(Stats.pctlerp(.25, f));
		vd.add(Stats.pctlerp(.50, f));
		vd.add(Stats.pctlerp(.75, f));
		break;
	    }
	}
	GlobalEntry g = new GlobalEntry(vs, vd, nat, mut, cwd, useCorresp);
	if(type == RunType.LOCAL)
	    outputSemiprocessedMutationLine(g, wt_local, outputMap);
	else
	    semiProcessed.add(g);
    }
    /* (non-Javadoc)
     * @see parse.TsParser#finish()
     */
    @Override
    public void finish() {
	// if GLOBAL:
	//   calculate quartiles over all wt lines
	//   for each stored quartile
	//     calculate and output quartile score w.r.t. wt
	if(type == RunType.GLOBAL) {
	    Vector<double[]> wt = processWT(wt_global);
	    for(GlobalEntry g : semiProcessed)
		outputSemiprocessedMutationLine(g, wt, outputMap);
	}
    }
    /**
     * Parses a series of lines split into fields and returns a vector of sorted arrays, one per numeric column.
     * @param wt_in input lines to process
     * @return sorted arrays, one per numeric column
     */
    protected Vector<double[]> processWT(Vector<String[]> wt_in) {
	Vector<double[]> wt = new Vector<double[]>();
	for(int i = 0; i < inputMap.length; i++) {
	    switch(inputMap[i]) {
	    case IGNORE:
	    case DONTCONVERT:
		break;
	    case CONVERT:
		// parse strings into array
		double f[] = new double[wt_in.size()];
		for(int j = 0; j < f.length; j++)
		    f[j] = Double.parseDouble(wt_in.get(j)[i]);
		Arrays.sort(f);
		wt.add(f);
		break;
	    }
	}
	return wt;
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
     * Does output on a partially processed set up mutation lines;
     * takes the place of outputMutationLine in V1 and V2 parsers.
     * @param g contains the partially processed mutation line elements
     * @param wt sorted list of wt entries for quartile calculations
     * @param map field map to use for output
     */
    public void outputSemiprocessedMutationLine(GlobalEntry g, Vector<double[]> wt, FieldType[] map) {
	String prot = Field.protP.matcher(g.stringElement.elementAt(0)).replaceAll("$1");
	int fi = 4;
	for(int i = 0; i < seqopt.length; i++) {
	    switch(seqopt[i]) {
	    case A:
		g.stringElement.add(fi, Field.chgFmt.format(Residue.getAminoChange(g.nativeRes, g.mutantRes)));
		break;
	    case B:
		g.stringElement.add(fi, Field.chgFmt.format(Residue.getAminoChangeGlenn(g.nativeRes, g.mutantRes)));
		break;
	    case C:
		g.stringElement.add(fi, Sequence.getAminoChangeID(g.nativeRes));
		break;
	    case D:
		g.stringElement.add(fi, Sequence.getAminoChangeID(g.mutantRes));
		break;
	    case E:
		g.stringElement.add(fi, Sequence.getAminoChangeGlennID(g.nativeRes));
		break;
	    case F:
		g.stringElement.add(fi, Sequence.getAminoChangeGlennID(g.mutantRes));
		break;
	    case G:
		g.stringElement.add(fi, Field.oddsFmt.format(Sequence.getAminoChangeOdds(g.nativeRes, g.mutantRes)));
		break;
	    case H:
		g.stringElement.add(fi, Field.oddsFmt.format(Sequence.getAminoChangeGlennOdds(g.nativeRes, g.mutantRes)));
		break;
	    case P:
		try {
		    g.stringElement.add(fi, Field.oddsFmt.format(seqan.getMutantLogOdds(g, prot, SUBDIR1)));
		} catch(SequenceException se) {
		    if(!badProt.contains(prot)) {
			System.err.println(se.getMessage());
			badProt.add(prot);
		    }
		    g.stringElement.add(fi, "?");
		}
		break;
	    case Q:
		try {
		    g.stringElement.add(fi, Field.oddsFmt.format(seqan.getNativeLogOdds(g, prot, SUBDIR1)));
		} catch(SequenceException se) {
		    if(!badProt.contains(prot)) {
			System.err.println(se.getMessage());
			badProt.add(prot);
		    }
		    g.stringElement.add(fi, "?");
		}
		break;
	    case R:
		try {
		    g.stringElement.add(fi, Field.oddsFmt.format(seqan.getNativeLogOdds(g, prot, SUBDIR1) - seqan.getMutantLogOdds(g, prot, SUBDIR1)));
		} catch(SequenceException se) {
		    if(!badProt.contains(prot)) {
			System.err.println(se.getMessage());
			badProt.add(prot);
		    }
		    g.stringElement.add(fi, "?");
		}
		break;
	    case S:
		try {
		    g.stringElement.add(fi, Field.chgFmt.format(seqan.getInfoContent(g, prot, SUBDIR1)));
		} catch(SequenceException se) {
		    if(!badProt.contains(prot)) {
			System.err.println(se.getMessage());
			badProt.add(prot);
		    }
		    g.stringElement.add(fi, "?");
		}
		break;
	    case T:
		try {
		    g.stringElement.add(fi, Field.defFmt.format(seqan.getMutantFrequency(g, prot, SUBDIR1)));
		} catch(SequenceException se) {
		    if(!badProt.contains(prot)) {
			System.err.println(se.getMessage());
			badProt.add(prot);
		    }
		    g.stringElement.add(fi, "?");
		}
		break;
	    case U:
		try {
		    g.stringElement.add(fi, Field.defFmt.format(seqan.getNativeFrequency(g, prot, SUBDIR1)));
		} catch(SequenceException se) {
		    if(!badProt.contains(prot)) {
			System.err.println(se.getMessage());
			badProt.add(prot);
		    }
		    g.stringElement.add(fi, "?");
		}
		break;
	    case V:
		try {
		    g.stringElement.add(fi, Field.defFmt.format(seqan.getNativeFrequency(g, prot, SUBDIR1) - seqan.getMutantFrequency(g, prot, SUBDIR1)));
		} catch(SequenceException se) {
		    if(!badProt.contains(prot)) {
			System.err.println(se.getMessage());
			badProt.add(prot);
		    }
		    g.stringElement.add(fi, "?");
		}
		break;
	    case p:
		try {
		    g.stringElement.add(fi, Field.oddsFmt.format(seqan.getMutantLogOdds(g, prot, SUBDIR2)));
		} catch(SequenceException se) {
		    if(!badProt.contains(prot)) {
			System.err.println(se.getMessage());
			badProt.add(prot);
		    }
		    g.stringElement.add(fi, "?");
		}
		break;
	    case q:
		try {
		    g.stringElement.add(fi, Field.oddsFmt.format(seqan.getNativeLogOdds(g, prot, SUBDIR2)));
		} catch(SequenceException se) {
		    if(!badProt.contains(prot)) {
			System.err.println(se.getMessage());
			badProt.add(prot);
		    }
		    g.stringElement.add(fi, "?");
		}
		break;
	    case r:
		try {
		    g.stringElement.add(fi, Field.oddsFmt.format(seqan.getNativeLogOdds(g, prot, SUBDIR2) - seqan.getMutantLogOdds(g, prot, SUBDIR2)));
		} catch(SequenceException se) {
		    if(!badProt.contains(prot)) {
			System.err.println(se.getMessage());
			badProt.add(prot);
		    }
		    g.stringElement.add(fi, "?");
		}
		break;
	    case s:
		try {
		    g.stringElement.add(fi, Field.chgFmt.format(seqan.getInfoContent(g, prot, SUBDIR2)));
		} catch(SequenceException se) {
		    if(!badProt.contains(prot)) {
			System.err.println(se.getMessage());
			badProt.add(prot);
		    }
		    g.stringElement.add(fi, "?");
		}
		break;
	    case t:
		try {
		    g.stringElement.add(fi, Field.defFmt.format(seqan.getMutantFrequency(g, prot, SUBDIR2)));
		} catch(SequenceException se) {
		    if(!badProt.contains(prot)) {
			System.err.println(se.getMessage());
			badProt.add(prot);
		    }
		    g.stringElement.add(fi, "?");
		}
		break;
	    case u:
		try {
		    g.stringElement.add(fi, Field.defFmt.format(seqan.getNativeFrequency(g, prot, SUBDIR2)));
		} catch(SequenceException se) {
		    if(!badProt.contains(prot)) {
			System.err.println(se.getMessage());
			badProt.add(prot);
		    }
		    g.stringElement.add(fi, "?");
		}
		break;
	    case v:
		try {
		    g.stringElement.add(fi, Field.defFmt.format(seqan.getNativeFrequency(g, prot, SUBDIR2) - seqan.getMutantFrequency(g, prot, SUBDIR2)));
		} catch(SequenceException se) {
		    if(!badProt.contains(prot)) {
			System.err.println(se.getMessage());
			badProt.add(prot);
		    }
		    g.stringElement.add(fi, "?");
		}
		break;
	    default:
		System.err.println("internal error: handling for SeqField " + seqopt[i] + " (" + seqopt[i].getFieldName() + ") not implemented");
		System.exit(1);
	    }
	    fi++;
	}
	int si = 0, qi = 0, wi = 0;  // stringElement, quartileElement, and wt indices, respectively
	for(int i = 0; i < map.length; i++) {
	    switch(map[i]) {
	    case IGNORE: break;
	    case DONTCONVERT:
		if(i > 0) System.out.format(",");
		System.out.format("%s", g.stringElement.get(si));
		si++;
		break;
	    case CONVERT:
		if(i > 0) System.out.format(",");
		System.out.print(Field.defFmt.format(Stats.normfind(g.quartileElement.get(qi), wt.get(wi)))
			+ "," + Field.defFmt.format(Stats.normfind(g.quartileElement.get(qi+1), wt.get(wi)))
			+ "," + Field.defFmt.format(Stats.normfind(g.quartileElement.get(qi+2), wt.get(wi))));
		qi += 3;
		wi++;
		break;
	    }
	}
	System.out.format("\n");
    }
    /**
     * Encapsulates information about a group of mutations lines that have been processed
     * to yield DONTCONVERT and CONVERT fields.  These are stored sequentially (without gaps for
     * IGNORE fields) in {@link GlobalEntry#stringElement} and {@link GlobalEntry#quartileElement},
     * respectively.
     * @author crispy
     */
    public static class GlobalEntry {
	public Vector<String> stringElement;
	public Vector<Double> quartileElement;
	public ResidueList nativeRes, mutantRes;
	public String cwd;
	public boolean useCorresp;
	public GlobalEntry(Vector<String> vs, Vector<Double> vd, ResidueList nat, ResidueList mut, String cwd, boolean corr) {
	    this.stringElement = vs;
	    this.quartileElement = vd;
	    this.nativeRes = nat;
	    this.mutantRes = mut;
	    this.cwd = cwd;
	    this.useCorresp = corr;
	}
    }
}
