/**
 * Created on Sep 1, 2009
 * Copyright Chris Poultney
 */
package parse;

import java.util.Vector;
import util.ResidueList;
import util.Field.FieldType;

/**
 * Encapsulates the implementation-independent aspects of producing an output file.
 * This will eventually be moved to a separate class/package.
 * @author crispy
 */
public interface TsParser {
    /**
     * Adds a parser-specific option.
     * @param args the command line args list
     * @param idx index of the arg to add
     * @return index of the next arg to parse
     * @throws ArgumentException on argument parsing error or usage request
     */
    public int registerOption(String[] args, int idx) throws ArgumentException;
    /**
     * Comma-separated list of column headers for extra fields to add (e.g. aminochange).
     * @return extra field names
     */
    public String getExtraFields();
    /**
     * Process and output the header line.
     * @param header header line, broken into fields
     * @param inputMap array of field types corresponding to header line fields
     */
    public void processHeader(String[] header, FieldType[] inputMap);
    /**
     * Any processing that needs to be done on a per-line (instead of per-residue) basis.
     * @param name current mutation description
     * @param mut current mutation(s)
     * @param fields the fields that make up the current line
     */
    public void processLine(String name, ResidueList mut, String[] fields);
    /**
     * Process a group of runs for wt or mutation.
     * @param name current mutation description (usually resi for single mutants)
     * @param mut AA mutation(s) (resn list)
     * @param nat native AA(s) (resn list)
     * @param currmut input lines; each entry is an array of fields, one entry per input line
     * @param isWt true for WT series, otherwise false
     * @param cwd is directory for original sample, used for loading additional data
     */
    public void processMutation(String name, ResidueList mut, ResidueList nat, Vector<String[]> currmut, boolean isWt, String cwd, boolean useCorresp);
    /**
     * Called when all input has been received.
     */
    public void finish();
}
