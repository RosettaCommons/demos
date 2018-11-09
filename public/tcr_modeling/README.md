TCRmodel Demo
=============

KEYWORDS: TCRmodel TCR modeling 

TCRmodel is a application that models T cell receptor (TCR) structures from sequence, developed by the Pierce Lab at the University of Maryland Institute for Bioscience and Biotechnology Research. Please see this reference for details on the server, and please cite it if you use the server in your work:

Gowthaman R, Pierce BG. (2018) TCRmodel: high resolution modeling of T cell receptors from sequence. Nucleic Acids Research, 46(Web Server issue), W396–W401. http://doi.org/10.1093/nar/gky432
Pubmed

Input files
-----------

* TCR sequences:

  TCR sequences for alpha and beta chain needs to be given as input separately. Use the flag -alpha and -beta to flags to provide the sequences for alpha and beta chains respectively. Here is an example sequence from the PDB (PDB:1TCR)
  Example alpha and beta chain sequence:
  -alpha QSVTQPDARVTVSEGASLQLRCKYSYSATPYLFWYVQYPRQGLQLLLKYYSGDPVVQGVNGFEAEFSKSNSSFHLRKASVHWSDSAVYFCAVSGFASALTFGSGTKVIVLPYIQNPEPAVYALKDPRSQDSTLCLFTDFDSQINVPKTMESGTFITDATVLDMKAMDSKSNGAIAWSNQTSFTCQDIFKETNATYPSSDVPC
  -beta EAAVTQSPRNKVAVTGGKVTLSCNQTNNHNNMYWYRQDTGHGLRLIHYSYGAGSTEKGDIPDGYKASRPSQENFSLILELATPSQTSVYFCASGGGGTLYFGAGTRLSVLEDLRNVTPPKVSLFEPSKAEIANKQKATLVCLARGFFPDHVELSWWVNGKEVHSGVSTDPQAYKESNYSYCLSSRLRVSATFWHNPRNHFRCQVQFHGLSEEDKWPEGSPKPVTQNISAEAWGRADC

* Template database:

  The template database should be placed by default in $Rosetta/database/additional_protocol_data/tcr. The template database can also be downloaded in to Rosetta/database using git clone. The command is: git clone git@github.com:RosettaCommons/additional_protocol_data.git
  The latest updated template database can be dowloaded from the TCrmodel web sever.  The link for the template database is: https://tcrmodel.ibbr.umd.edu/links
  If the template database is located in outside $Rosetta/database/additional_protocol_data/, then use the flag "-tcr_template_db_path" to give the path.

  				   				 
Running TCRmodel
----------------
  The TCR model application can be run by simply providing the alpha and beta chain sequences. Example command to run the TCR modeling application:

```
    $>  $ROSETTA3/bin/tcr.default.linuxgccrelease @tcr.options  
```

Refinement of the model:						       
    To refine the grafted model we can do the fullatom minimization separately in Rosetta or as an
additional option while running the modeling job.
  -minimize_model true/false [Default:true]
  -relax_model true/false [Default:false]

```
    $>  $ROSETTA3/bin/tcr.default.linuxgccrelease @tcr_relax.options  
```

  To run TCR modeling with optional user provided templates we need to specify the template files for the corresponding TCR segment. It is possible to provide the pdb structure file as template for any TCR segment. Also we can just provide the PDB ID as template, th ecorresponding PDB structure should be present in the template database. 

The TCR segments for this modeling purpose is grouped in to either one of the following groups. 
  # Framework, CDR1, CDR2 and CDR3 
  # Germline, CDR3

  Example flags for user provided template PDB id's
  -alpha_germline_template_id 1mwa_A
  -alpha_cdr3_template_id 2oi9_B
  -alpha_orientation_template_id 1mwa_A
  -beta_germline_template_id 2q86_B
  -beta_cdr3_template_id 2apb_A
  -beta_orientation_template_id 1mwa_B

Modeling with user provided PDB template id's:
```
    $>  $ROSETTA3/bin/tcr.default.linuxgccrelease @tcr_template_id.options  
```

  Example flags for user provided template PDB structures
  -alpha_germline_template_pdb gma_tmplt_piece.pdb
  -alpha_cdr3_template_pdb cdr3a_tmplt_piece.pdb
  -alpha_orientation_template_pdb ora_tmplt_piece.pdb
  -beta_germline_template_pdb gmb_tmplt_piece.pdb
  -beta_cdr3_template_pdb cdr3b_tmplt_piece.pdb
  -beta_orientation_template_pdb orb_tmplt_piece.pdb 

Modeling with user provided template structures:
```
    $>  $ROSETTA3/bin/tcr.default.linuxgccrelease @tcr_template_pdb.options  
```

Loop modeling and refinement:
  Loop modeling and refinement can be performed on the initial template base model structure. Loop modeling can be carried out for single CDR3 loop of alpha or beta chain or for both chains. 
  Example flags for the TCR modeling with loop modeling and refinement:
  -refine_tcr_cdr3_loops true/false [Refine the CDR3 loops of Alpha and Beta chain.Default:False]
  -remodel_tcr_cdr3_loops true/false [Remodel the CDR3 loops of Alpha and Beta chain.Default:False]
  -remodel_tcr_cdr3a_loop true/false [Remodel the CDR3 loop of alpha chain. Useful if remodeling is required only for the CDR3 loop of alpha chain. Default:False]
  -remodel_tcr_cdr3b_loop true/false [Remodel the CDR3 loop of beta chain. Useful if remodeling is required only for the CDR3 loop of beta chain.Default:False]
  -refine_all_tcr_cdr_loops true/false [Refine all the CDR loops. Refinement includes CDR1, CDR2, CDR3 & HV4 loops of Alpha and Beta chains. Default:False]

```
    $>  $ROSETTA3/bin/tcr.default.linuxgccrelease @tcr_loop.options  
```
 	
Analyzing the results
---------------------

  Expected output files are:
  Running TCR modeling application will result in the following PDB format output files. The -out::prefix flag can be used to add prefix to the output file names.
    tcrmodel.pdb The final refined final model.
    tcr_graftmodel.pdb The crude grafted model without refinement.
    tcr_loopmodel.pdb The model after application of optional loop modeling or refinement.

  The output log should look similar to:

  ~/Rosetta/main/source/bin/tcr.macosclangrelease -alpha QSVTQPDARVTVSEGASLQLRCKYSYSATPYLFWYVQYPRQGLQLLLKYYSGDPVVQGVNGFEAEFSKSNSSFHLRKASVHWSDSAVYFCAVSGFASALTFGSGTKVIVLPYIQNPEPAVYALKDPRSQDSTLCLFTDFDSQINVPKTMESGTFITDATVLDMKAMDSKSNGAIAWSNQTSFTCQDIFKETNATYPSSDVPC -beta EAAVTQSPRNKVAVTGGKVTLSCNQTNNHNNMYWYRQDTGHGLRLIHYSYGAGSTEKGDIPDGYKASRPSQENFSLILELATPSQTSVYFCASGGGGTLYFGAGTRLSVLEDLRNVTPPKVSLFEPSKAEIANKQKATLVCLARGFFPDHVELSWWVNGKEVHSGVSTDPQAYKESNYSYCLSSRLRVSATFWHNPRNHFRCQVQFHGLSEEDKWPEGSPKPVTQNISAEAWGRADC
  –––––––––––––––––––––––
  Output_model:					 tcrmodel.pdb		Score	-151.325
  Tcr Alpha truncated Domain sequence: VTQPDARVTVSEGASLQLRCKYSYSATPYLFWYVQYPRQGLQLLLKYYSGDPVVQGVNGFEAEFSKSNSSFHLRKASVHWSDSAVYFCAVSGFASALTFGSGTKVIVL
  Tcr Alpha Framework sequence : VTQPDARVTVSEGASLQLRCWYVQYPRQGLQLLLRKASVHWSDSAVYFCFGSGTKVIVL
  Tcr Alpha Domain CDR1 sequence: KYSYSATPYLF
  Tcr Alpha Domain CDR2 sequence: LKYYSGDPVVQGVNGFEAEFSKSNSSFH
  Tcr Alpha Domain CDR3 sequence: AVSGFASALT
  Tcr Alpha germline sequence : VTQPDARVTVSEGASLQLRCKYSYSATPYLFWYVQYPRQGLQLLLKYYSGDPVVQGVNGFEAEFSKSNSSFHLRKASVHWSDSAVYFCFGSGTKVIVL
  Alpha germline template : 1mwa_A_A
  Alpha CDR3 template : 2oi9_B_A
  Tcr Beta truncated Domain : VTQSPRNKVAVTGGKVTLSCNQTNNHNNMYWYRQDTGHGLRLIHYSYGAGSTEKGDIPDGYKASRPSQENFSLILELATPSQTSVYFCASGGGGTLYFGAGTRLSVL
  Tcr Beta Framework sequence : VTQSPRNKVAVTGGKVTLSCWYRQDTGHGLRLILILELATPSQTSVYFCFGAGTRLSVL
  Tcr Beta Domain CDR1 sequence: NQTNNHNNMY
  Tcr Beta Domain CDR2 sequence: HYSYGAGSTEKGDIPDGYKASRPSQENFS
  Tcr Beta Domain CDR3 sequence: ASGGGGTLY
  Tcr Beta germline sequence : VTQSPRNKVAVTGGKVTLSCNQTNNHNNMYWYRQDTGHGLRLIHYSYGAGSTEKGDIPDGYKASRPSQENFSLILELATPSQTSVYFCFGAGTRLSVL
  Beta germline template : 2q86_B_B
  Beta CDR3 template : 2apb_A_B
  Alpha Beta orientation template : 1mwa_A_A 1mwa_B_B

  In this case the Rosetta score for the final modeled pose is -151.325.

Other options for TCRmodel:

* '-template_similarity_cutoff' : Similarity cutoff to ignore similar template sequences from template database. Default:100
* '-template_identity_cutoff' : Identity cutoff to ignore similar template sequences from template database.Default:100
* '-blastp_identity_cutoff' : Identity cutoff to ignore similar template sequences from template database.Default:100
* '-ignore_list' : List of PDB id's to ignore as templates.Default:None Default=False 	Boolean
* '-use_beta_germline_templates' : Use germline templates for beta chain, by default germline or framework templates choosen by sequence match.
