/**
 * Created on Sep 1, 2009
 * Copyright Chris Poultney
 */
package test;

import java.util.Arrays;
import util.Stats;

/**
 * Should give same results as ml/tmp.d/percentile-check.ods.
 * @author crispy
 */
public class QuartileTest2 {
    public static void main(String[] args) {
	double wt[] =  { 0.916, 0.926, 0.926, 0.931, 0.931, 0.931, 0.936, 0.936, 0.936, 0.936,
			 0.941, 0.941, 0.941, 0.946, 0.946, 0.946, 0.946, 0.946, 0.946, 0.946,
			 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950,
			 0.950, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955,
			 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.965, 0.970 };
	double mut[] = { 0.931, 0.946, 0.950, 0.936, 0.950, 0.941, 0.931, 0.936, 0.916, 0.926,
			 0.941, 0.926, 0.941, 0.936, 0.931, 0.931, 0.921, 0.950, 0.931, 0.916,
			 0.941, 0.926, 0.960, 0.901, 0.931, 0.960, 0.941, 0.936, 0.950, 0.946,
			 0.926, 0.946, 0.921, 0.921, 0.960, 0.931, 0.926, 0.936, 0.941, 0.936,
			 0.941, 0.941, 0.946, 0.936, 0.960, 0.946, 0.926, 0.911, 0.941, 0.926 };
	Arrays.sort(wt);
	Arrays.sort(mut);
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
	    System.out.format("v%.6f = %.6f, %.6f\n", x, p, Stats.normfind(p, wt));
	}
    }
}
