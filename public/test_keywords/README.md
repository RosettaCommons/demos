Keyword Testing (Internal detail)
=================================

This directory is our test to make sure that every tutorial/demo has a keyword list and that the keyword is in the approved list.
Each tutorial and demo should have keyword list as follows:

	KEYWORDS: STRUCTURE_PREDICTION LOOPS GENERAL DESIGN MEMBRANES NUCLEIC_ACIDS RNA DNA PEPTIDES COMPARATIVE_MODELING  

Keywords one and two denote organizational levels.  If an item does not have a level2, use general.  Any extra keywords are used for keyword search.
The approved keywords can be found in keywords.txt in the demo root directory.  Keywords are present in order to have one approved name for something, such as SYMMETRY vs SYMMETRIC or DOCK vs DOCKING.  

The following command will test every directory in this type, via integration.py on the test server.

	$> python $ROSETTA_DEMOS/check_demos_for_keywords.py ../ $ROSETTA_DEMOS/keywords.txt

