/**
 * Created on Sep 16, 2009
 * Copyright Chris Poultney
 */
package util;


import java.util.Arrays;
import java.util.Random;

/**
 * Single class for all statistics methods.
 * @author crispy
 */
public class Stats {
    /**
     * Returns the value at the given percentile, interpolating between values if necessary.
     * @param p the percentile to find, range [0-1]
     * @param a the SORTED array to search
     * @return the interpolated value
     */
    public static double pctlerp(double p, double[] a) {
	// from wikipedia: for p=0..100, n=1..N: n=(p/100)(N-1)+1
	// for me: p=0..1, n=0..N-1: n=p(N-1)
	int l = a.length - 1;
	if(p <= 0) return a[0];
	else if(p >= 1) return a[l];
	else {
	    double n = p * l;		// floating point position in array
	    int i = (int)Math.floor(n);	// integer component of n
	    double d = n - i;		// fractional component of n
	    return a[i] + d * (a[i+1] - a[i]);
	}
    }
    /**
     * Finds the index of v within a, normalized and clamped to the range 0-1.
     * @param v the value to search for
     * @param a the SORTED array to search
     * @return the normalized, clamped index value
     */
    public static double normfind(double v, double[] a) {
	int l = a.length - 1;
	if(v <= a[0]) return 0;
	else if(v >= a[l]) return 1;
	else {
	    int i = Arrays.binarySearch(a, v);
	    if(i > 0) {		// found match; for consistency, find earliest match
		while(i > 0 && a[i-1] == a[i])
		    i--;
		return i / (double)l;
	    } else {		// no match found; must interpolate between nearest values
		int ip = -i-2;	// closest array index before v: -(i+1)-1
		double n = ip / (double)l;	// percentile of lower value
		return n + ((v-a[ip]) / (a[ip+1]-a[ip])) / (double)l;
	    }
	}
    }
    /**
     * Returns statistics on a list of numbers
     * @param f input array
     * @return array [ mean, median, variance ]
     */
    public static double[] getStats(double[] f, boolean median_compat) {
	double v[] = new double[3];
	Arrays.sort(f);
	// median
	if(f.length % 2 == 0 && median_compat)
	    v[1] = (f[f.length/2-1] + f[f.length/2])/2.0;
	else
	    v[1] = f[f.length/2];
	// mean & variance
	double m = 0, ms = 0;
	for(int i = 0; i < f.length; i++) {
	    m += f[i];
	    ms += f[i] * f[i];
	}
	v[0] = m / (double)f.length;  // mean
	v[2] = ms / (double)f.length - v[0]*v[0];  // variance
	return v;
    }
    /**
     * Calculate how many times a process with chance of success <prob_random> succeeds more than
     * <num_actual> times in <num_pred> independent samples; perform this process ntrials times,
     * and report the success fraction as p-value.
     * @param ntrials
     * @param prob_rnd
     * @param num_pred
     * @param num_act
     * @return
     */
    public static double pValue(int ntrials, double prob_rnd, int num_pred, int num_act) {
	Random rnd = new Random();
	int success = 0;
	for(int i = 0; i < ntrials; i++) {
	    int n = 0;
	    for(int j = 0; j < num_pred; j++) {
		if(rnd.nextDouble() < prob_rnd)
		    n++;
	    }
	    if(n > num_act)
		success++;
	}
	return (double)success / ntrials;
    }
    /**
     * Calculates the fraction of trials for which a random process outperforms the process which
     * produced the given results.  For each trial, num_pred samples are drawn, without replacement,
     * from an initial pool of cnt_all samples.  cnt_des of these samples have the desired label,
     * implicitly giving an initial probability cnt_des/cnt_all of drawing a desired sample at random.
     * A trial outperforms the original process if num_act or more of the num_pred samples drawn have
     * the desired label.
     * @param ntrials
     * @param cnt_all
     * @param cnt_des
     * @param num_pred
     * @param num_act
     * @return
     */
    public static double pValueNoRep(int ntrials, int cnt_all, int cnt_des, int num_pred, int num_act) {
	Random rnd = new Random();
	int success = 0;
	int[] pool = new int[cnt_all];
	{
	    int i;
	    for(i = 0; i < cnt_des; i++)
		pool[i] = 1;
	    for( ; i < cnt_all; i++)
		pool[i] = 0;
	}
	// a little tricky to speed things up:
	//  sz is the number of samples left in the pool
	//  when a sample is drawn, it is swapped with the last sample in the pool (index sz-1)
	//  sz is then decremented
	//  this avoids having to reinitialize the array at the start of each trial
	for(int i = 0; i < ntrials; i++) {
	    int sz = pool.length;
	    int n = 0;
	    for(int j = 0; j < num_pred; j++) {
		int x = rnd.nextInt(sz);
		int v = pool[x];
		if(v == 1)
		    n++;
		pool[x] = pool[sz-1];
		pool[sz-1] = v;
		sz--;
	    }
	    if(n >= num_act)
		success++;
	}
	return (double)success / ntrials;
    }
}
