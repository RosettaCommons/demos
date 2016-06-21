package util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.Hashtable;
import java.util.TreeSet;
import parse.ParseV2;
import util.Residue.AA;

public class Sequence {
    static final boolean debug = false;
    static final String AC_FN = "data/ac.txt";
    static int[][] acOdds;
    static final String ACG_FN = "data/acg.txt";
    static int[][] acgOdds;
    static {
	try {
	    parseAminoChangeTable();
	} catch(Exception e) {
	    System.err.println("could not parse aminochange table " + AC_FN + ", error follows");
	    e.printStackTrace();
	    System.exit(1);
	}
	try {
	    parseAminoChangeGlennTable();
	} catch(Exception e) {
	    System.err.println("could not parse aminochangeglenn table " + ACG_FN + ", error follows");
	    e.printStackTrace();
	    System.exit(1);
	}
    }
    public static enum SeqField {
	A("aminochange"), B("aminochangeglenn"), C("acfrom"), D("acto"), E("acgfrom"), F("acgto"),
	G("acodds"), H("acgodds"),
	P("pssm_mut"), Q("pssm_nat"), R("pssm_diff"), S("info_cont"),
	T("freq_mut"), U("freq_nat"), V("freq_diff"),
	p("pssm_mut2"), q("pssm_nat2"), r("pssm_diff2"), s("info_cont2"),
	t("freq_mut2"), u("freq_nat2"), v("freq_diff2");
	private final String name;
	SeqField(String n) { this.name = n; }
	public String getFieldName() { return name; }
    }
    public static SeqField getField(char c) { return getField(c+""); }
    public static SeqField getField(String s) {
	return SeqField.valueOf(s);
    }
    public static String getAminoChangeID(ResidueList r) {
	if(r.aa.length != 1)
	    return "?";
	else
	    return Residue.getAAProp(r.aa[0]).id() + "";
    }
    public static String getAminoChangeGlennID(ResidueList r) {
	if(r.aa.length != 1)
	    return "?";
	else
	    return Residue.getAAPropGlenn(r.aa[0]).id() + "";
    }
    public static int getAminoChangeOdds(AA nat, AA mut) {
	int ni = Residue.getAAProp(nat).id(), mi = Residue.getAAProp(mut).id();
	return acOdds[ni][mi];
    }
    public static double getAminoChangeOdds(ResidueList nat, ResidueList mut) {
	double o = 0;
	for(int i = 0; i < nat.aa.length; i++)
	    o += getAminoChangeOdds(nat.aa[i], mut.aa[i]);
	return o / nat.aa.length;
    }
    public static int getAminoChangeGlennOdds(AA nat, AA mut) {
	int ni = Residue.getAAPropGlenn(nat).id(), mi = Residue.getAAPropGlenn(mut).id();
	return acgOdds[ni][mi];
    }
    public static double getAminoChangeGlennOdds(ResidueList nat, ResidueList mut) {
	double o = 0;
	for(int i = 0; i < nat.aa.length; i++)
	    o += getAminoChangeGlennOdds(nat.aa[i], mut.aa[i]);
	return o / nat.aa.length;
    }
    public static class SequenceAnalyzer {
	Hashtable<String,PSSM> pssmH = new Hashtable<String,PSSM>();
	TreeSet<String> pssmBad = new TreeSet<String>();
	Hashtable<String,AACorresp> correspH = new Hashtable<String,AACorresp>();
	TreeSet<String> correspBad = new TreeSet<String>();
	String MODELLER_ROOT;
	public SequenceAnalyzer() {
	    MODELLER_ROOT = System.getProperty("modeller.root");
	    if(MODELLER_ROOT ==  null) {
		MODELLER_ROOT = "none";
		System.err.println("Could not determine MODELLER_ROOT, many sequence-based terms will be unavailable");
		System.err.println("  to enable, launch with \"-Dmodeller.root=$MODELLER_ROOT\"");
	    }
	}
	public AACorresp getAACorresp(ParseV2.GlobalEntry g, String prot) throws IOException {
	    // throw exception on first error, return null on subsequent errors
	    if(MODELLER_ROOT.equals("none"))
		return null;
	    // get corresp file, or add to "bad" list and bail
	    AACorresp cc;
	    if(!g.useCorresp)
		cc = new AACorresp();
	    else {
		String path = MODELLER_ROOT + "/" + g.cwd + "/" + prot + ".crs";
		try {
		    cc = AACorresp.makeCorresp(path);
		} catch(IOException ie) {
		    correspBad.add(path);
		    throw new IOException("CRS file error: " + path, ie);
		}
	    }
	    return cc;
	}
	public PSSM getPSSM(ParseV2.GlobalEntry g, String prot, String subdir) throws IOException {
	    if(MODELLER_ROOT.equals("none"))
		return null;
	    // get pssm, or add to "bad" list and bail
	    String path = MODELLER_ROOT + "/" + g.cwd + "/" + subdir + "/" + prot + ".cp.gz";
	    PSSM pssm = pssmH.get(path);
	    if(pssm == null) {
		if(!pssmBad.contains(path))
		    try {
//			System.err.println("loading " + path);
			pssm = PSSM.parsePSSM(new URL("file:" + path));
			pssmH.put(path, pssm);
//			System.err.println("  done");
		    } catch(IOException ie) {
			pssmBad.add(path);
			ie.printStackTrace();
			throw new IOException("PSSM file error, PSSM terms will not be available: " + path, ie);
		    }
	    }
	    return pssm;
	}
	public double getNativeLogOdds(ParseV2.GlobalEntry g, String prot, String subdir) throws SequenceException {
	    try {
		AACorresp cc = getAACorresp(g, prot);
		PSSM pssm = getPSSM(g, prot, subdir);
		if(cc == null || pssm == null) {
		    throw new SequenceException("CRS or PSSM unavailable: " + prot);
		} else {
		    double sc = 0;
		    for(int a = 0; a < g.nativeRes.aa.length; a++) {
			if(g.nativeRes.pos == null)
			    throw new Exception("position data not available");
			int seqposNat = cc.translate(g.nativeRes.pos[a]);
			if(!g.nativeRes.aa[a].getShort().name().equals(pssm.getSequenceAt(seqposNat)+""))
			    throw new Exception("NLO: native res does not match pssm seq: " + prot + ", mutation " + g.mutantRes);
			sc += pssm.getScore(seqposNat, g.nativeRes.aa[a].getShort().name());
		    }
		    sc /= (double)g.nativeRes.aa.length;
		    return sc;
		}
	    } catch(Exception e) {
		throw new SequenceException(e);
	    }
	}
	public double getMutantLogOdds(ParseV2.GlobalEntry g, String prot, String subdir) throws SequenceException {
	    try {
		AACorresp cc = getAACorresp(g, prot);
		PSSM pssm = getPSSM(g, prot, subdir);
		if(cc == null || pssm == null) {
		    throw new SequenceException("CRS or PSSM unavailable: " + prot);
		} else {
		    double sc = 0;
		    for(int a = 0; a < g.mutantRes.aa.length; a++) {
			if(g.mutantRes.pos == null)
			    throw new Exception("position data not available");
			int seqposNat = cc.translate(g.nativeRes.pos[a]), seqposMut = cc.translate(g.mutantRes.pos[a]);
			if(!g.nativeRes.aa[a].getShort().name().equals(pssm.getSequenceAt(seqposNat)+""))
			    throw new Exception("MLO: native res does not match pssm seq: " + prot + " " +
				    g.nativeRes.aa[a].getShort().name() + g.nativeRes.pos[a] + " != " + pssm.getSequenceAt(seqposNat) +
				    " (a=" + a + ", seqposNat=" + seqposNat + ", pos=" + g.nativeRes.pos[a] + ")");
			sc += pssm.getScore(seqposMut, g.mutantRes.aa[a].getShort().name());
		    }
		    sc /= (double)g.mutantRes.aa.length;
		    return sc;
		}
	    } catch(Exception e) {
		throw new SequenceException(e);
	    }
	}
	public double getInfoContent(ParseV2.GlobalEntry g, String prot, String subdir) throws SequenceException {
	    try {
		AACorresp cc = getAACorresp(g, prot);
		PSSM pssm = getPSSM(g, prot, subdir);
		if(cc == null || pssm == null) {
		    throw new SequenceException("CRS or PSSM unavailable: " + prot);
		} else {
		    double sc = 0;
		    for(int a = 0; a < g.nativeRes.aa.length; a++) {
			if(g.nativeRes.pos == null)
			    throw new Exception("position data not available");
			int seqposNat = cc.translate(g.nativeRes.pos[a]);
			if(!g.nativeRes.aa[a].getShort().name().equals(pssm.getSequenceAt(seqposNat)+""))
			    throw new Exception("IC: native res does not match pssm seq: " + prot + ", mutation " + g.mutantRes);
			sc += pssm.getInfoContent(seqposNat);
		    }
		    sc /= (double)g.nativeRes.aa.length;
		    return sc;
		}
	    } catch(Exception e) {
		throw new SequenceException(e);
	    }
	}
	public double getNativeFrequency(ParseV2.GlobalEntry g, String prot, String subdir) throws SequenceException {
	    try {
		AACorresp cc = getAACorresp(g, prot);
		PSSM pssm = getPSSM(g, prot, subdir);
		if(cc == null || pssm == null) {
		    throw new SequenceException("CRS or PSSM unavailable: " + prot);
		} else {
		    double sc = 0;
		    for(int a = 0; a < g.nativeRes.aa.length; a++) {
			if(g.nativeRes.pos == null)
			    throw new Exception("position data not available");
			int seqposNat = cc.translate(g.nativeRes.pos[a]);
			if(!g.nativeRes.aa[a].getShort().name().equals(pssm.getSequenceAt(seqposNat)+""))
			    throw new Exception("NF: native res does not match pssm seq: " + prot + ", mutation " + g.mutantRes);
			sc += pssm.getWeightedFreq(seqposNat, g.nativeRes.aa[a].getShort().name());
		    }
		    sc /= (double)g.nativeRes.aa.length;
		    return sc;
		}
	    } catch(Exception e) {
		throw new SequenceException(e);
	    }
	}
	public double getMutantFrequency(ParseV2.GlobalEntry g, String prot, String subdir) throws SequenceException {
	    try {
		AACorresp cc = getAACorresp(g, prot);
		PSSM pssm = getPSSM(g, prot, subdir);
		if(cc == null || pssm == null) {
		    throw new SequenceException("CRS or PSSM unavailable: " + prot);
		} else {
		    double sc = 0;
		    for(int a = 0; a < g.mutantRes.aa.length; a++) {
			if(g.mutantRes.pos == null)
			    throw new Exception("position data not available");
			int seqposNat = cc.translate(g.nativeRes.pos[a]), seqposMut = cc.translate(g.mutantRes.pos[a]);
			if(!g.nativeRes.aa[a].getShort().name().equals(pssm.getSequenceAt(seqposNat)+""))
			    throw new Exception("MF: native res does not match pssm seq: " + prot + ", mutation " + g.mutantRes);
			sc += pssm.getWeightedFreq(seqposMut, g.mutantRes.aa[a].getShort().name());
		    }
		    sc /= (double)g.mutantRes.aa.length;
		    return sc;
		}
	    } catch(Exception e) {
		throw new SequenceException(e);
	    }
	}
    }
    public static void parseAminoChangeTable() throws IOException, NumberFormatException, ClassNotFoundException {
	acOdds = new int[4][4];
	URL u = getResource(AC_FN);
	BufferedReader in = new BufferedReader(new InputStreamReader(u.openStream()));
	String s;
	int l = 0;
	while((s = in.readLine()) != null) {
	    String[] f = s.split("\\s+");
	    for(int i = 0; i < 4; i++)
		acOdds[l][i] = Integer.parseInt(f[i]);
	    l++;
	}
	in.close();
	if(debug) {
	    System.err.println("aminochange table");
	    for(int m = 0; m < acOdds.length; m++) {
		for(int n = 0; n < acOdds[m].length; n++)
		    System.err.format("%d\t", acOdds[m][n]);
		System.err.println();
	    }
	}
    }
    public static void parseAminoChangeGlennTable() throws IOException, NumberFormatException, ClassNotFoundException {
	acgOdds = new int[7][7];
	URL u = getResource(ACG_FN);
	BufferedReader in = new BufferedReader(new InputStreamReader(u.openStream()));
	String s;
	int l = 0;
	while((s = in.readLine()) != null) {
	    String[] f = s.split("\\s+");
	    for(int i = 0; i < 7; i++)
		acgOdds[l][i] = Integer.parseInt(f[i]);
	    l++;
	}
	in.close();
	if(debug) {
	    System.err.println("aminochangeglenn table");
	    for(int m = 0; m < acgOdds.length; m++) {
		for(int n = 0; n < acgOdds[m].length; n++)
		    System.err.format("%d\t", acgOdds[m][n]);
		System.err.println();
	    }
	}
    }
    public static URL getResource(String s) throws ClassNotFoundException {
	ClassLoader l = util.Sequence.class.getClassLoader();
        return l.getResource(s);
    }
    public static class SequenceException extends Exception {
	private static final long serialVersionUID = 1L;
	public SequenceException(Throwable t) { super(t); }
	public SequenceException(String m) { super(m); }
	public SequenceException(String m, Throwable t) { super(m, t); }
    }
}
