/**
 * Created on Sep 1, 2009
 * Copyright Chris Poultney
 */
package util;

import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.regex.Pattern;

/**
 * @author crispy
 *
 */
public class Field {
    /** derive protein name from full protein+mutation string */
    public static final Pattern protP = Pattern.compile("(^.*)-[^-]*$");
    public static final String[] IGNORE_FIELDNAMES = { "ara", "description", "id",
	"fa_atr", "fa_sol", "fa_intra_rep", "irms", "maxsub", "rms" };
    public static final String[] DONTCONVERT_FIELDNAMES = { "file", "RES", "amino", "MUT", "SPECIES", "aminochange", "aminochangeglenn",
	"acfrom", "acto", "acgfrom", "acgto", "acodds", "acgodds",
	"pssm_mut", "pssm_nat", "pssm_diff", "info_cont", "freq_mut", "freq_nat", "freq_diff",
	"pssm_mut2", "pssm_nat2", "pssm_diff2", "info_cont2", "freq_mut2", "freq_nat2", "freq_diff2",
//	"ref", "dslf_ss_dih",
	"PT", "ss_S", "ss_H", "ss_L", "ACC", "ACCF", "ACCFSA", "ACCP", "ACCPSA", "ACCPR", "ACCPRSA", "PT_val" };
    public static enum FieldType { CONVERT, IGNORE, DONTCONVERT };
    public static DecimalFormat defFmt = new DecimalFormat("#####.###");
    public static DecimalFormat chgFmt = new DecimalFormat("#.##");
    public static DecimalFormat oddsFmt = new DecimalFormat("#.#");
    /**
     * Assigns field handling instructions according to {@link #IGNORE_FIELDNAMES}, etc.
     * @param header array of text header fields
     * @return mapping of field number (0-indexed) to handling type
     */
    public static FieldType[] getFieldMap(String[] header) {
	Arrays.sort(IGNORE_FIELDNAMES);
	Arrays.sort(DONTCONVERT_FIELDNAMES);
	FieldType[] m = new FieldType[header.length];
	for(int i = 0; i < m.length; i++)
	    if(Arrays.binarySearch(IGNORE_FIELDNAMES, header[i]) >= 0)
		m[i] = FieldType.IGNORE;
	    else if(Arrays.binarySearch(DONTCONVERT_FIELDNAMES, header[i]) >= 0)
		m[i] = FieldType.DONTCONVERT;
	    else
		m[i] = FieldType.CONVERT;
	return m;
    }
}
