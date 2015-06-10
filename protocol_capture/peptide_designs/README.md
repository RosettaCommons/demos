### General Information ##################
# Your name: Deanne Sammond

# Protocol Name: Peptide backbone/sequence design


# Brief Description:
This protocol builds (or extends) a backbone for a peptide bound to a target protein, then designs a low-energy sequence. 

# if this is a protocol attached to a rosettaCon talk provide:

  # Title, Authors & Lab , Year, Session and Day of talk
  Computational design of a new protein-protein interface between Gi1 and a redesigned RGS14 GoLoco, Deanne Sammond, Dustin Bosch, Glenn Butterfoss, Mischa Machius, David Siderovski, Brian Kuhlman, Kuhlman lab, 2010, Session 3 on Wednesday August 4th
  # Abstract
Our project is the computational design of a new high-affinity protein-protein interface.  Our model system is an x-ray crystal structure of  Gi1 bound to the GoLoco domain from the RGS14 protein.  RGS14 GoLoco spans two domains of Gi1, with the C-terminal random coil region binding to the all-helical domain of Gi1.  We removed this C-terminal portion of GoLoco, replacing the random coil with a de novo designed alpha helix.  The redesigned GoLoco binds to Gi1 with a dissociation constant of 810nM, the correct binding of the newly designed GoLoco was confirmed using disruptive mutations at the Gi1:GoLoco interface, and the correctness of the computational design was assessed with by x-ray crystallography.

### running #########
# Important Flags, other flags that can be used, required flag, brief description of each flag or option:
-ex1, -ex2, -exOH, -extrachi_cutoff 1 all seem to be very important for the sequence design run.

# Example Rosetta Command Line:
This protocol uses 2 separate rosetta runs - one is centroid mode to build backbone coordinates, and the other is a design run to find a low-energy sequence
 rosetta.mactel aa input_pdb _ -s g000.pdb -loops

 rosetta.mactel -design -l list_of_pdbs -tail -begin 342 -end 351 -chain_ -series bb -protein g000 -resfile g000_resfile -ex1 -ex2 -extrachi_cutoff 1 -exOH -no_his_his_pairE -tight_hb -try_both_his_tautomers 

# Example Overall Command Line (if overall protocol is run via a script or other program)
The entire protocol uses 2 calls to Rosetta and 2 additional scripts.  For a complete description of how I ran the protocol see README_dws_in_detail.

### versions #########
# If checked into trunk, svn revision number: 
29304, https://svn.rosettacommons.org/source/trunk/rosetta++
# If NOT checked into trunk, path and svn revision number: 

# Version for other codes used, version for interfaces or libs used (if relevant)


# References to published works using this protocol (you may include submitted, and in press works)
Manuscript in preparation.

# Other Comments: 



