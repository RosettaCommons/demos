/**
 * Created on Oct 18, 2009
 * Copyright Chris Poultney
 */
package app;

import util.Stats;

/**
 * @author crispy
 * Simple class to calculate p-values.
 */
public class PVal {
    public static void usage() { usage(null); }
    public static void usage(String msg) {
	if(msg != null)
	    System.err.println(msg);
	System.err.println("usage: java app.PVal <ntrials> <prob_random> <num_pred> <num_actual>");
	System.err.println("Calculates what fraction of times a random process outperforms prediction results");
	System.err.println("  ntrials: how many groups of num_pred samples to test");
	System.err.println("  prob_random: the chance of success of an independent random observation");
	System.err.println("  num_pred: how many samples to draw in each trial");
	System.err.println("  num_actual: number of successful samples from prediction");
	System.exit(1);
    }
    /**
     * @param args
     * Calculates p-value; see {@link #usage(String)} for details.
     */
    public static void main(String[] args) {
	if(args.length != 4)
	    usage();
	// get args
	int ntrials = getInt(args[0], "ntrials");
	double prob_rnd = getDouble(args[1], "prob_random");
	int num_pred = getInt(args[2], "num_pred");
	int num_act = getInt(args[3], "num_actual");
	// calculate p-value
	double p = Stats.pValue(ntrials, prob_rnd, num_pred, num_act);
	System.out.format("p-value: %.3g\n", p);
    }
    public static int getInt(String s, String desc) {
	try {
	    int n = Integer.parseInt(s);
	    return n;
	} catch(Exception e) {
	    usage("error: could not parse value of argument " + desc);
	    return 0;
	}
    }
    public static double getDouble(String s, String desc) {
	try {
	    double x = Double.parseDouble(s);
	    return x;
	} catch(Exception e) {
	    usage("error: could not parse value of argument " + desc);
	    return 0;
	}
    }
}
