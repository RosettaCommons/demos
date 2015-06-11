/**
 * Created on Aug 12, 2009
 * Copyright Chris Poultney
 */
package app;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.util.Vector;
import java.util.regex.Pattern;
import parse.ArgumentException;
import parse.ParseV1;
import parse.ParseV2Old;
import parse.ParseV2;
import parse.TsParser;
import util.Field;
import util.ResidueList;
import util.Field.FieldType;

/**
 * Takes a *-in.csv file and converts it to Weka input csv, replacing parsedectree.
 * @author crispy
 */
public class TsParseCSV {
    String infile;
    TsParser parser;
    /**
     * Records file name and options for {@link TsParseCSV#convert()}.
     * @param infile input CSV file
     */
    public TsParseCSV(String infile, TsParser parser) {
	this.infile = infile;
	this.parser = parser;
    }
    /**
     * Main function loop: converts input CSV to output CSV.  Input assumptions: contains header line;
     * first three columns are name, mutation, and label.
     * @throws IOException
     */
    public void convert() throws IOException {
	LineNumberReader in;
	if(infile.equals("-"))
	    in = new LineNumberReader(new BufferedReader(new InputStreamReader(System.in)));
	else
	    in = new LineNumberReader(new BufferedReader(new FileReader(infile)));
	String line;
	int lineno = 0;
	Pattern subP = Pattern.compile("^#sub,(.*)$");
	Pattern addP = Pattern.compile("^([^,]+,[^,]+,[^,]+,[^,]+)");
	String lastName = "", lastPt = "",  lastProt = "", cwd = "", lastcwd = "";
	boolean useCorresp = false, lastCorresp = false;
	ResidueList lastMut = null, nat = null, lastNat = null;
	Vector<String[]> currmut = new Vector<String[]>();  // each inner String[] is a series of fields from currline below
	while((line = in.readLine()) != null) {
	    lineno++;
	    // header line
	    if(line.startsWith("RES")) {
		line = subP.matcher(line).replaceAll("$1");
		String[] f = line.split(",");
		FieldType inputMap[] = Field.getFieldMap(f);
		line = addP.matcher(line).replaceAll("$1," + parser.getExtraFields());
		f = line.split(",");
		parser.processHeader(f, inputMap);
		continue;
	    } else if(line.startsWith("#")) {
		if(line.startsWith("##dir:")) {
		    cwd = line.substring(6);
		    useCorresp = false;
		} else if(line.startsWith("##corresp"))
		    useCorresp = true;
		continue;
	    }
	    // all other lines
	    String[] f = line.split(",");
	    String name = f[0], pt = f[2];  // these first three are assumed
	    ResidueList mut = new ResidueList(f[1]);
	    String prot = Field.protP.matcher(f[0]).replaceAll("$1");
	    parser.processLine(name, mut, f);
	    if(pt.equals("NAT"))   // for now, ignore NAT lines
		continue;
	    if(!lastProt.equals(prot)) {
		;
	    }
	    if(!lastName.equals(name)) {
		// end-of-mutation-block (name) processing
		// start new block (name)
		if(!pt.equals("wt")) {
		    System.err.format("First group for new name must be wt (received \"%s\"), line %s%n", pt, in.getLineNumber());
		    usage();
		}
		nat = mut;
	    }
	    if(lastMut != null && (!lastProt.equals(prot) || !lastName.equals(name) || !lastMut.equals(mut))) {
		// end-of-mutation processing
		parser.processMutation(lastName, lastMut, lastNat, currmut, lastPt.equals("wt"), lastcwd, lastCorresp);
		if(currmut.size() != 50)
		    throw new IOException("unexpected number of input lines: " + currmut.size() + ", line " + in.getLineNumber());
		// start new mutation
		currmut.clear();
	    }
	    currmut.add(f);
	    lastName = name;
	    lastMut = mut;
	    lastPt = pt;
	    lastProt = prot;
	    lastNat = nat;
	    lastcwd = cwd;
	    lastCorresp = useCorresp;
	}
	in.close();
	parser.processMutation(lastName, lastMut, lastNat, currmut, lastPt.equals("wt"), lastcwd, useCorresp);
	parser.finish();
    }
    /**
     * Prints usage message and exits; unified exit point for any exception or error condition.
     */
    public static void usage() { usage(""); }
    /**
     * Prints usage message plus parser-specific info; unified exit point for any exception or error condition.
     * @param versionSpecific parser-specific usage information
     */
    public static void usage(String versionSpecific) {
	System.err.println("usage: TsParseCSV <-v1|-v2|-v2old> [parser-options] <file-in.csv>");
	System.err.println("converts *-in.csv file to Weka input csv file");
	System.err.println("  <file-in.csv> can be \"-\" to specify stdin");
	System.err.print(versionSpecific);
	System.exit(1);
    }
    /**
     * Launches TsParseCSV; type java TsParseCSV --usage for command line args.
     * @param args
     */
    public static void main(String[] args) {
	TsParser p = null;
	try {
	    int i = 0;
	    while(i < args.length && args[i].startsWith("-")) {
		if(args[i].equals("-v1")) {
		    p = new ParseV1();
		    i++;
		} else if(args[i].equals("-v2old")) {
		    p = new ParseV2Old();
		    i++;
		} else if(args[i].equals("-v2")) {
		    p = new ParseV2();
		    i++;
		} else if(args[i].equals("-")) {	// read from standard input; should be last argument
		    break;
		} else {
		    if(p == null) {
			System.err.format("unrecognized global option %s, aborting\n", args[0]);
			usage();
		    } else {
			i = p.registerOption(args, i);
		    }
		}
	    }
	    if(p == null) {
		System.err.println("no parser specified, use -v<n>");
		usage();
	    }
	    if(args.length == i) {
		System.err.println("input file not provided");
		usage();
	    }
	    TsParseCSV parse = new TsParseCSV(args[i], p);
	    parse.convert();
	} catch(ArgumentException ae) {
	    usage(ae.getMessage());
	} catch(Exception e) {
	    e.printStackTrace();
	    usage();
	}
    }
}
