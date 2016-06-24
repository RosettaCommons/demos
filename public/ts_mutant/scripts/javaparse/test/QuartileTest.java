/**
 * Created on Sep 1, 2009
 * Copyright Chris Poultney
 */
package test;

import util.Stats;

/**
 * Quick test of quartile methods.
 * @author crispy
 */
public class QuartileTest {
    public static void main(String[] args) {
	double wt[] = { 54, 58, 60, 62, 67, 70, 73, 74, 75, 78, 79, 80, 81, 85, 85, 86, 87, 87, 88, 88 };
	double mut[] = { 73, 76, 77, 77, 78, 80, 84, 87, 87, 89 };
	System.out.println("wt0: " + Stats.pctlerp(.0, wt));
	System.out.println("wt25: " + Stats.pctlerp(.25, wt));
	System.out.println("wt50: " + Stats.pctlerp(.50, wt));
	System.out.println("wt75: " + Stats.pctlerp(.75, wt));
	System.out.println("wt100: " + Stats.pctlerp(1.0, wt));
	System.out.println("wt/wt");
	eval(wt, wt);
	System.out.println("mut/wt");
	eval(mut, wt);
    }
    public static void eval(double[] mut, double[] wt) {
	double v[] = { 0.0, 0.25, 0.5, 0.75, 1.0 };
	for(double x: v) {
	    double p = Stats.pctlerp(x, mut);
	    System.out.format("v%.2f = %.2f, %.4f\n", x, p, Stats.normfind(p, wt));
	}
    }
}
