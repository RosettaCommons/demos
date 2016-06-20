rna_denovo.macosgccrelease -database  ~/minirosetta_database/ -fasta gcuucggc.fasta -nstruct 1 -out:file:silent gcuucggc.out -minimize_rna -cycles 5000 -mute all -native gcuucggc_RNA.pdb  > farfar.log


rna_denovo.macosgccrelease -database  ~/minirosetta_database/ -fasta gcuucggc.fasta -nstruct 1 -out:file:silent gcuucggc_NATIVE.out -minimize_rna -cycles 5000 -mute all -native gcuucggc_RNA.pdb  -vall_torsions 1f7y.torsions > farfar_native.log

