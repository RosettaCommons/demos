package app;

import util.Stats;

/**
 * Variant of app.PVal that draws draws samples without replacement from an starting pool.
 * @author crispy
 */
public class PVal2 {
    public static void usage() { usage(null); }
    public static void usage(String msg) {
	if(msg != null)
	    System.err.println(msg);
	System.err.println("usage: java app.PVal2 <ntrials> <cnt_all_samples> <cnt_desired_samples> <num_pred> <num_actual>");
	System.err.println("Calculates what fraction of times a random process outperforms prediction results");
	System.err.println("  ntrials: how many groups of num_pred samples to test");
	System.err.println("  cnt_all_samples: total number of samples");
	System.err.println("  cnt_desired_samples: number of samples with the desired label");
	System.err.println("  num_pred: number of samples predicted to have the desired label (total predicted)");
	System.err.println("  num_actual: number of successful samples from prediction (correct predicted)");
	System.exit(1);
    }
    /**
     * @param args
     */
    public static void main(String[] args) {
	if(args.length != 5)
	    usage();
	// get args
	int ntrials = PVal.getInt(args[0], "ntrials");
	int cnt_all = PVal.getInt(args[1], "cnt_all_samples");
	int cnt_des = PVal.getInt(args[2], "cnt_desired_samples");
	int num_pred = PVal.getInt(args[3], "num_pred");
	int num_act = PVal.getInt(args[4], "num_actual");
	// calculate p-value
	double p = Stats.pValueNoRep(ntrials, cnt_all, cnt_des, num_pred, num_act);
	if(p == 0)
	    System.out.format("p-value: <%.3g\n", (1.0/ntrials));
	else
	    System.out.format("p-value: %.3g\n", p);
    }
}
