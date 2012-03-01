awk '   BEGIN { start1 = 0; start2 = 0; cnt = 0;
		ONE="AcCDEFGHIKLMNPQRSTVWY";
		THREE="ALA.cys.CYS.ASP.GLU.PHE.GLY.HIS.ILE.LYS.LEU.MET.ASN.PRO.GLN.ARG.SER.THR.VAL.TRP.TYR"
		print ">";
	}
     	NF~1 {  if($1 == "_Mol_residue_sequence") start1 = 1;
		else if($1 == ";") {
			if(start1 == 1 && start2 == 0) start2 = 1;
			else if(start1 == 1 && start2 == 1)
			{ start1 = 0; start2 = 0;}
		}
		else if( start2==1) {
			printf("%s", $1);
		}
	}
	END { print; }
	' $1

# CS-ROSETTA: System for Chemical Shifts based protein structure prediction using ROSETTA
# (C) Shen and Bax 2007-2008, Lab of Chemical Physics, NIDDK, NIH
# Version 1.01(build 2009.1117.15)
#
# bmrb2fasta.com: search BMRB header for fasta sequence
#         syntax: bmrb2fasta.com bmrb.str > seq.fasta
