# AbPredict

KEYWORDS: STRUCTURE_PREDICTION ANTIBODIES

**Summary**: Methods for antibody structure prediction rely on sequence homology to experimentally determined structures. Resulting models may be accurate, but they are often stereochemically strained, limiting their usefulness in modeling and design workflows. We present the AbPredict 2 webserver, which instead of using sequence homology, conducts a Monte Carlo based search for low-energy combinations of backbone conformations to yield accurate and unstrained antibody structures.  
**Availability & implementation**: We introduce several important improvements over the previous AbPredict implementation: (1) backbones and sidechains are now modeled using ideal bond lengths and angles, reducing stereochemical strain to levels observed in high-resolution experimental structures; (2) sampling of the rigid-body orientation at the light-heavy chain interface is improved, increasing model accuracy; and (3) runtime is reduced 20-fold without compromising accuracy, enabling the implementation of AbPredict 2 as a fully automated webserver (http://abpredict.weizmann.ac.il). Accurate and unstrained antibody model structures may in some cases obviate the need for experimental structures in antibody optimization workflows.

For additional reading see: [27717001](https://www.ncbi.nlm.nih.gov/pubmed/27717001)

The AbPredict XML (AbPredict_xsd.xml) uses Rosetta's Splice movers to generate a model of an Fv. 

The input sequence should be passed as a script var in the following syntax, either as a command line argument or added to the flags file:

    -parser:script_vars sequence=IKMTQSPSSMYASLGERVTITCKASQDIRKYLNWYQQKPWKSPKTLIYYATSLADGVPSRFSGSGSGQDYSLTISSLESDDTATYYCLQHGESPYTFGGGTKLEIQLQQSGAELVRPGALVKLSCKASGFNIKDYYMHWVKQRPEQGLEWIGLIDPENGNTIYDPKFQGKASITADTSSNTAYLQLSSLTSEDTAVYYCARDNSYYFDYWGQGTTLTVS 

**Note:** that the sequence should only be that of the Fv fragment and formatted in the following way:
* The Vl sequence comes before the Vh sequence
2. The N-terminus tail length (the residues before the first disulfide cysteine, not including) of the light chain should be exactly 21 aa's long
3. The The N-terminus tail length of the heavy chain (starting from the first residue after the Vl/Vh chain break up to the first disulfide cys, not including) should be exactly 20 aa's
4. There should be exactly 7 aa's after the conserved L3 phe (L98 in Chothia numbering)
5. There should be exactly 0 aa's after the conserved H3 trp (H103 in Chothia numbering)


In order to run _AbPredict_ we need to select fragments of the appropriate length from the conformation database. 
This can be done by running the following bash script:
`./create_run.sh <VL length> <L3 length> <HL length> <L3 length>`

*   VL length/HL length : The number of aa's between the disulfide cysteines (including)
*   L3 length/H3 length : The number of aa's from the second disulfide cysteine (including) to the conserved L3/H3 phe/trp.


The script will generate a file named "segment_lengths_script_vars", each line in the file containts script vars for a single simulation trajectory. e.g.:


    -parser:script_vars entry_H1_H2=3OAUH entry_L1_L2=1H8OA entry_H3=4NHGD entry_L3=1NGXL

To run a single simulation trajectory:

    <Rosetta_Directory>/main/source/bin/rosetta_scripts.default.linuxgccrelease @flags -parser:script_vars entry_H1_H2=2IBZX entry_L1_L2=3DSFL entry_H3=3V4UH entry_L3=1LO2L  sequence=IKMTQSPSSMYASLGERVTITCKASQDIRKYLNWYQQKPWKSPKTLIYYATSLADGVPSRFSGSGSGQDYSLTISSLESDDTATYYCLQHGESPYTFGGGTKLEIQLQQSGAELVRPGALVKLSCKASGFNIKDYYMHWVKQRPEQGLEWIGLIDPENGNTIYDPKFQGKASITADTSSNTAYLQLSSLTSEDTAVYYCARDNSYYFDYWGQGTTLTVS 

The output pdb will be saved in the pdb directory. It is recommended to run a few hundred trajectories per single modeling task. The final models should be clustered and the best Rosetta energy structure from the top 3 clusters chosen for continued processing. 



