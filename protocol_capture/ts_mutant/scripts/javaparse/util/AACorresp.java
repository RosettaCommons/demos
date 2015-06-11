package util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.Map;
import java.util.TreeMap;

/**
 * Handles conversion between (possibly discontinuous) residue numbers in PDB file
 * and continuous 1-indexed numbering used by PSI-BLAST, etc.
 * @author crispy
 */
public class AACorresp {
    TreeMap<Integer,Integer> corresp;
    public AACorresp() {
	corresp = new TreeMap<Integer,Integer>();
	corresp.put(new Integer(1), new Integer(1));
    }
    private AACorresp(TreeMap<Integer,Integer> corresp) {
	this.corresp = corresp;
    }
    public int translate(int from) {
	Map.Entry<Integer,Integer> m = corresp.floorEntry(new Integer(from));
	int offset = m.getKey().intValue() - m.getValue().intValue();
	return from - offset;
    }
    public static AACorresp makeCorresp(String path) throws IOException {
	TreeMap<Integer,Integer> cr = new TreeMap<Integer,Integer>();
	LineNumberReader in = new LineNumberReader(new BufferedReader(new FileReader(path)));
	String s;
	while((s = in.readLine()) != null) {
	    if(s.startsWith("seq")) {
		;  // header line
	    } else {
		try {
		    String[] f = s.split("\\s+");
		    int from = Integer.parseInt(f[1]), to = Integer.parseInt(f[0]);
		    cr.put(new Integer(from), new Integer(to));
		} catch(Exception e) {
		    throw new IOException("error parsing .crs file " + path + " at line " + in.getLineNumber());
		}
	    }
	}
	return new AACorresp(cr);
    }
}
