#!/usr/bin/tcsh

set PATH_TO_DEMO=/vol/ek/ravehb/rosetta/RosettaCon2011/demos/global_dock_ssrA_peptide_against_sspB/
set PATH_TO_SCRIPTS=$PATH_TO_DEMO/scripts

cd $PATH_TO_DEMO/output_files
mkdir Clustering
cd Clustering
$PATH_TO_SCRIPTS/clustering/cluster.sh 500 2 ../score.sc ../native.pdb ../decoys.silent reweighted_sc
$PATH_TO_SCRIPTS/clustering/join_clustering_info.sh > join.txt
head -n 11 join.txt

