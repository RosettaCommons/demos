homology_modeling_with_end_extension README
Created: Aug 5, 2011
Last Updated: Aug 5, 2011

This demo "homology_modeling_with_end_extension" will illustrate a protocol for modeling the structure of a protein given its sequence, provided:
- there are sequence homologs in the PDB, but
- an n-terminal region of around 140 residues has no identifiable sequence similarity to a protein structure in the pdb

The input to the protocol is a sequence of interest.  The output will be a set of Rosetta-generated structural models.

Outline
We will break down homology_modeling_with_end_extension into 2 parts:
Part I describes a protocol for homology modeling of the protein (the region of the protein that has a homolog in the pdb).
Part II describes a protocol for ab initio modeling of the n-terminal region while keeping the previously modeled homologous region in Part 1 rigid.

Part I - homology modeling

There are 2 ways to do the homology modeling:
Option 1) use the rosetta_comparative_modeling_protocol or
Option 2) use the Robetta server (http://www.robetta.org).  Robetta can automatically generate models for the region that has a homolog in the pdb. If using Robetta, use its output (the generated models for the homologous region) and skip to Part II (Broker Rigid Chunk Claimer).

The following 7 steps are for the rosetta_comparative_modeling_protocol option, based on a protocol developed by James Thompson and TJ Brunette.

1. get homologs and alignments.

Prior to homology modeling, PBD homologs for the protein of interest need to be identified and alignments generated.  Several external programs can be used for this step like the HHpred server (http://toolkit.tuebingen.mpg.de/hhpred), psi-blast, or hhsearch.

In this demo, we will use HHpred. To submit an HHpred job go to http://toolkit.tuebingen.mpg.de/hhpred, paste the sequence of your protein of interest into the text field, and then click the "Submit job" button. Getting results should take a few minutes. On the results page, you can save the alignments to a file by clicking the Save link under the "Results" tab. We include an example output file below:

<pre>
./starting_files/query.hhr
</pre>

At this point, you may truncate the query so that the long (not well aligned) n-terminal region is removed. Then resubmit the truncated fasta to HHpred for more accurate alignments and continue with the comparative modeling protocol using the truncated sequence. For this example, we do not truncate the sequence for simplicity.  We will model the N-terminal region in Part II using the the broker in context of the homology model.


2. get template PDBs.

run
<pre>
./scripts/cm_scripts/bin/get_pdbs.pl ./starting_files/query.hhr
</pre>

This script will get the PDB for each alignment from the RCSB and clean the PDB files so they are compatible with Rosetta. The PDB files are named by the PDB code and chain id followed by a ".pdb" file extension.


3. convert alignment format to 'grishin' (Rosetta uses this format.)

run
<pre>
./scripts/cm_scripts/bin/convert_aln.pl -format_in hhsearch -format_out grishin ./starting_files/query.hhr > query.grishin
</pre>

This script will generate the following output file:
query.grishin

The 'grishin' format is the standard alignment format that is used by Rosetta. An example of this format is shown below:

<pre>
## t288_ 1be9A_4
# hhsearch_3 33
scores_from_program: 0.000000 0.998400
2 VPGKVTLQKDAQNLIGISIGGGAQYCPCLYIVQVFDNTPAALDGTVAAGDEITGVNGRSIKGKTKVEVAKMIQEVKGEVTIHYNKLQ
9 EPRRIVIHRGS-TGLGFNIIGGED-GEGIFISFILAGGPADLSGELRKGDQILSVNGVDLRNASHEQAAIALKNAGQTVTIIAQYKP
--
</pre>

The first two lines represent the identifier of the query and template sequences, The string "1be9A_4" should be a unique identifier for the sequence alignment where the first five characters are the PDB code and chain id of the template PDB. The starting positions are given at the begining of each alignment line.


4. fix alignments to match the PDBs (HHpred alignments may not match the sequence and/or amino acid positions in the PDBs).  This is common and should be fixed using the following application. The example command below assumes that the Rosetta release is installed in the directory where this demo exists.

run
<pre>
../rosetta-3.3/rosetta_source/bin/fix_alignment_to_match_pdb.macosgccrelease -in:file:alignment query.grishin -cm:aln_format grishin -out:file:alignment query_fixed.grishin -database ../rosetta-3.3/rosetta_database -in:file:template_pdb *.pdb
</pre>

This script will generate the following output file:
query_fixed.grishin


5. generate fragment files

The next step is to generate 3mer and 9mer fragment files, which can easily be done using the Robetta fragment server at http://robetta.bakerlab.org/fragmentsubmit.jsp.

Save the fragment files along with the psipred secondary structure prediction file locally.

We include fragment files and the psipred secondary structure prediction file:

<pre>
./rosetta_inputs/aat000_03_05.200_v1_3  (3mer fragment file)
./rosetta_inputs/aat000_09_05.200_v1_3  (9mer fragment file)
./rosetta_inputs/t000_.psipred_ss2      (psipred file)
</pre>

6. run the rosetta comparative modeling application (see documentation for the comparative modeling application for more information)

The Rosetta comparative modeling application 1) generates an incomplete model based on the template structure by copying coordinates over the aligned regions, 2) rebuilds the missing parts using the Rosetta loop modeling protocol, and then 3) runs full-atom refinement of the protein model using the Rosetta full-atom energy function. For more information, see the documentation for the "loop modeling" and "relax" applications.

The example command below uses the following input files:

<pre>
./starting_files/query.fasta                 (query fasta sequence file)
./rosetta_inputs/aat000_03_05.200_v1_3       (3mer fragment file)
./rosetta_inputs/aat000_09_05.200_v1_3       (9mer fragment file)
./rosetta_inputs/t000_.psipred_ss2           (psipred file)
./query_fixed.grishin                        (the alignment file)
*.pdb                                        (the template PDB files)
</pre>

<pre>
../rosetta-3.3/rosetta_source/bin/minirosetta.macosgccrelease @./rosetta_inputs/comparative_modeling.args -loops:frag_files ./rosetta_inputs/aat000_09_05.200_v1_3 ./rosetta_inputs/aat000_03_05.200_v1_3 none -in:file:psipred_ss2 ./rosetta_inputs/t000_.psipred_ss2 -in:file:alignment ./rosetta_inputs/query_fixed.grishin -in:file:template_pdb *.pdb -database ../rosetta-3.3/rosetta_database -out:file:silent query.out -in:file:fasta ./starting_files/query.fasta -nstruct 1
</pre>

This application will generate the following output file which contains the models in a compressed "silent file" format:
query.out

For simplicity we set -nstruct to 1 in this demo, but in real world applications you will want to set it higher to incease the amount of conformational sampling. The standard protocol is to generate around 1,000 to 10,000 separate models, select the lowest 10% of models by Rosetta energy, and then choose clusters using the "Cluster" application.


7. Extract pdbs from silent file:

<pre>
../rosetta-3.3/rosetta_source/bin/extract_pdbs.macosgccrelease -in:file:fullatom -in:file:silent query.out -database ../rosetta-3.3/rosetta_database
</pre>

This script will generate extacted PDB files with the following name format where 1N9LA would be replaced with the PDB code and chain id of the template:

S_1N9LA_1_0001.pdb


Select a homology model to use in Part II. In real world applications you would base selection on the results of clustering as briefly described above.


Part II: Broker Rigid Chunk Claimer (Bruno Correia's protocol)

Given the homology model, Part II describes how to model the terminal region using ab initio modeling in the context of the homology model generated in Part I.  To do this, we will use the fragment insertion sampling protocol with the "broker rigid chunk claimer" in Rosetta.

1. Generate a full length model if the n-terminal region was omitted

run
<pre>
../rosetta-3.3/rosetta_source/bin/full_length_model.macgccrelease -in:file:fasta ./starting_files/query.seq -loops:frag_files ./rosetta_inputs/aat000_09_05.200_v1_3 ./rosetta_inputs/aat000_03_05.200_v1_3 none -loops:frag_sizes 9 3 1 -in:file:s S_1N9LA_1_0001.pdb
</pre>

2. Make a region file that specifies the region you'd like to keep rigid (this is the region that was modeling in Part I using the homology modeling application). The format of the region file consists of a single line that contains the start and end positions that span the region that will be kept rigid.  In the example below positions 143-236 will be kept rigid.

example region file:
./rosetta_inputs/1N9LA_143_236.region

format of the regions file:
<pre>
RIGID 143 236 0 0 0
</pre>

3. Create a setup file for using the RigidChunkClaimer

example setup file:
./rosetta_inputs/1N9LA.tpb

format of the setup file:
<pre>
CLAIMER RigidChunkClaimer
REGION_FILE 1N9LA_143_236.region
PDB_FILE 1N9LA.pdb
END_CLAIMER
</pre>

4. Run the broker protocol using Rosetta:

../rosetta-3.3/rosetta_source/bin/minirosetta.macosgccrelease @./rosetta_inputs/broker.args -in:file:fasta ./starting_files/query.fasta -broker:setup ./rosetta_inputs/1N9LA.tpb -database ../rosetta-3.3/rosetta_database -frag3 ./rosetta_inputs/aat000_03_05.200_v1_3 -frag9 ./rosetta_inputs/aat000_09_05.200_v1_3 -nstruct 1

For the simplicity of this protocol demo, we set -nstruct to 1. In real world applications, we suggest making many more, perhaps 10000 model structures, dependent on the length of the extension.

5. Evaluate the models via clustering.  Refer to the clustering documentation for more information.





