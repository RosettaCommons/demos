/**
 * Created on Aug 13, 2009
 * Copyright Chris Poultney
 */
package util;

import java.util.Arrays;
import java.util.List;

/**
 * @author crispy
 *
 */
public class Residue {
    public static enum AAShort {
	A(AA.ALA), C(AA.CYS), D(AA.ASP), E(AA.GLU), F(AA.PHE), G(AA.GLY), H(AA.HIS), I(AA.ILE), K(AA.LYS), L(AA.LEU),
	M(AA.MET), N(AA.ASN), P(AA.PRO), Q(AA.GLN), R(AA.ARG), S(AA.SER), T(AA.THR), V(AA.VAL), W(AA.TRP), Y(AA.TYR);
	private final AA al;
	AAShort(AA a) { al = a; }
	public AA getLong() { return al; }
    }
    public static enum AA {
	ALA, CYS, ASP, GLU, PHE, GLY, HIS, ILE, LYS, LEU,
	MET, ASN, PRO, GLN, ARG, SER, THR, VAL, TRP, TYR;
	public AAShort getShort() {
	    for(AAShort a : AAShort.values())
		if(a.getLong() == this)
		    return a;
	    return null;
	}
    }
    public static enum AAProp {
	SmallHydrophobic(0, new AA[] { AA.GLY, AA.ALA }),
	LargeHydrophobic(1, new AA[] { AA.VAL, AA.LEU, AA.ILE, AA.MET, AA.PRO, AA.PHE, AA.TRP }),
	Polar(2, new AA[] { AA.SER, AA.THR, AA.ASN, AA.GLN, AA.CYS, AA.TYR }),
	Charged(3, new AA[] { AA.LYS, AA.ARG, AA.HIS, AA.ASP, AA.GLU });
	private final List<AA> res;
	private final int _id;
	AAProp(int _id, AA[] rl) { this._id = _id; res = Arrays.asList(rl); }
	public List<AA> getList() { return res; }
	public int id() { return _id; }
    };
    public static enum AAPropGlenn {
	G1(0, new AA[] { AA.GLY, AA.PRO }),
	G2(1, new AA[] { AA.VAL, AA.ALA, AA.ILE, AA.LEU }),
	G3(2, new AA[] { AA.MET, AA.PHE, AA.TRP, AA.TYR }),
	G4(3, new AA[] { AA.LYS, AA.ARG, AA.HIS }),          // + charge
	G5(4, new AA[] { AA.GLN, AA.ASN }),
	G6(5, new AA[] { AA.ASP, AA.GLU }),                  // - charge
	G7(6, new AA[] { AA.SER, AA.THR, AA.CYS });          // small polar
	private final List<AA> res;
	private final int _id;
	AAPropGlenn(int _id, AA[] rl) { this._id = _id; res = Arrays.asList(rl); }
	public List<AA> getList() { return res; }
	public int id() { return _id; }
    }
    public static AAProp getAAProp(AA r) {
	for(AAProp prop : AAProp.values())
	    if(prop.getList().contains(r))
		return prop;
	return null;
    }
    public static AAPropGlenn getAAPropGlenn(AA r) {
	for(AAPropGlenn prop : AAPropGlenn.values())
	    if(prop.getList().contains(r))
		return prop;
	return null;
    }
    public static int getAminoChange(AA nat, AA mut) {
	if(nat == mut)
	    return 0;
	else if(getAAProp(nat) == getAAProp(mut))
	    return 1;
	else
	    return 2;
    }
    public static int getAminoChangeGlenn(AA nat, AA mut) {
	if(nat == mut)
	    return 0;
	else if(getAAPropGlenn(nat) == getAAPropGlenn(mut))
	    return 1;
	else
	    return 2;
    }
    public static double getAminoChange(ResidueList a, ResidueList b) {
	double x = 0;
	for(int i = 0; i < a.aa.length; i++)
	    x += getAminoChange(a.aa[i], b.aa[i]);
	return x / (double)a.aa.length;
    }
    public static double getAminoChangeGlenn(ResidueList a, ResidueList b) {
	double x = 0;
	for(int i = 0; i < a.aa.length; i++)
	    x += getAminoChangeGlenn(a.aa[i], b.aa[i]);
	return x / (double)a.aa.length;
    }
    public static void check() {
	for(AA a : AA.values())
	    if(a.getShort().getLong() != a)
		System.out.println(a + ": match error");
	    else
		System.out.println(a + ": OK");
	for(AAShort a : AAShort.values()) {
	    if(a.getLong().getShort() != a)
		System.out.println(a + ": match error");
	    else
		System.out.println(a + ": OK");
	}
	System.out.println();
	System.out.format("%4s", "");
	for(AA b : AA.values())
	    System.out.format("%4s", b);
	System.out.println();
	for(AA a : AA.values()) {
	    System.out.format("%4s", a);
	    for(AA b : AA.values())
		System.out.format("%4d", getAminoChange(a, b));
	    System.out.println();
	}
	System.out.println();
	System.out.format("%4s", "");
	for(AA b : AA.values())
	    System.out.format("%4s", b);
	System.out.println();
	for(AA a : AA.values()) {
	    System.out.format("%4s", a);
	    for(AA b : AA.values())
		System.out.format("%4d", getAminoChangeGlenn(a, b));
	    System.out.println();
	}
    }
}
