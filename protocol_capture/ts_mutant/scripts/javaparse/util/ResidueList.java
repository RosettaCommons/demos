/**
 * Created on Sep 15, 2009
 * Copyright Chris Poultney
 */
package util;

import util.Residue.AA;

/**
 * Centralized implementation of a list of residues for handling multiple-residue mutants.
 * @author crispy
 */
public class ResidueList {
    public AA[] aa;
    public int[] pos;
    public ResidueList(String s) {
	String[] m = s.split("\\|");
	String[] f = m[0].split(":");
	aa = new AA[f.length];
	for(int i = 0; i < f.length; i++)
	    aa[i] = AA.valueOf(f[i]);
	pos = null;
	if(m.length > 1) {
	    f = m[1].split(":");
	    pos = new int[f.length];
	    for(int i = 0; i < f.length; i++)
		pos[i] = Integer.parseInt(f[i]);
	}
    }
    public boolean equals(Object other) {
	if(this == other) return true;
	else if(other == null) return false;
	else if(this.getClass() != other.getClass()) return false;
	ResidueList otherList = (ResidueList)other;
	if(this.aa.length != otherList.aa.length)
	    return false;
	else
	    for(int i = 0; i < aa.length; i++)
		if(this.aa[i] != otherList.aa[i])
		    return false;
	return true;
    }
    public String toString() {
	if(aa.length == 0) {
	    return "(none)";
	} else {
	    StringBuffer s = new StringBuffer();
	    s.append(aa[0].name() + pos[0]);
	    for(int i = 1; i < aa.length; i++)
		s.append(":" + aa[i].name() + pos[i]);
	    return s.toString();
	}
    }
}
