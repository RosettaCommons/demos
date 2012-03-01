awk '   BEGIN { start0 = 0; start1 = 0; start2 = 0; cnt = 0;
		ONE   = "AcCDEFGHIKLMNPQRSTVWY";
		THREE = "ALA.cys.CYS.ASP.GLU.PHE.GLY.HIS.ILE.LYS.LEU.MET.ASN.PRO.GLN.ARG.SER.THR.VAL.TRP.TYR"
	}
	NF~9 && ($9 == "1"||$9 == "2"||$9 == "."){
		pos = index( THREE,toupper($4) );
		if( pos>0 && ($5=="CA"||$5=="CB"||$5=="C"||$5=="N"||$5=="H"||$5=="HN"||index($5,"HA")>0) ){
			atomName = $5;
			if( atomName == "H" ) atomName = "HN"
			printf(" %4d %1s %4s %8.3f\n",$3,substr(ONE,pos/4+1,1),atomName,$7);
		}
	}
	NF~8 && ($8 == "1"||$8 == "2"||$8 == "."){
		pos = index( THREE,toupper($3) );
		if( pos>0 && ($4=="CA"||$4=="CB"||$4=="C"||$4=="N"||$4=="H"||$4=="HN"||index($4,"HA")>0) ){
			atomName = $4;
			if( atomName == "H" ) atomName = "HN"
			printf(" %4d %1s %4s %8.3f\n",$2,substr(ONE,pos/4+1,1),atomName,$6);
		}
	}
     	NF~1 {  
		if($1 == "_Mol_residue_sequence") {
			start1 = 1;
		}
		else if($1 == ";") { 
			if(start1 == 1 && start2 == 0) start2 = 1;
			else if(start1 == 1 && start2 == 1) { 
				start1 = 0; start2 = 0; 
				if( start0 == 0) {
					print "\n\nVARS   RESID RESNAME ATOMNAME SHIFT "; 
					print "FORMAT %4d %1s %4s %8.3f\n";
					start0 = 1;
				}
			}
		}
		else if( start2 == 1 && start0 == 0 ) { 
			if( cnt == 0 ) printf("\nDATA SEQUENCE ");
			printf("%s ", substr($1,0,10));
			cnt = cnt + 1;
			if( cnt >= 5 ) cnt = 0;
			if( cnt == 0 ) printf("\nDATA SEQUENCE ");
			printf("%s ", substr($1,11,10));
			cnt = cnt + 1;
			if( cnt >= 5 ) cnt = 0;
		}
	}' $1

# CS-ROSETTA: System for Chemical Shifts based protein structure prediction using ROSETTA
# (C) Shen and Bax 2007-2008, Lab of Chemical Physics, NIDDK, NIH
#
# bmrb2talos.com: converts BMRB format chemical shift file to talos/mfr/csrosetta format
#         syntax: bmrb2talos.com bmrb.str > talos.tab
