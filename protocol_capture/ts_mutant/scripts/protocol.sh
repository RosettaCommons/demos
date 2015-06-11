# this makes the next steps easier, it’s the PDB file name without the .pdb
PN=YBR109C
# generate accessibility file
gen-acc ${PN}.pdb > sasa-${PN}.txt
# generate a master list of mutations to make using accessibility cutoff of 0
# cutoff of 10 is standard, we use 0 to reduce the number of runs is this example
gen_csv.sh -protein $PN -species Scer -cutoff 0
# make pdb file for each mutation
mkmut-csv ${PN}.csv
# generate relax run script for each mutation and native
# define $MINI_BIN and $MINI_DB as the Rosetta binary and database directories
#   on the machine where you will be performing the Rosetta runs
gen.sh -bin $MINI_BIN -db $MINI_DB -native ${PN}.pdb ${PN}-*.pdb

# do runs: this is system-dependent, as it usually requires access to a cluster
# the line below is just an example
for a in *.sh; do qsub -d $(pwd) $a; done

# add extra score terms from relax runs to final score output
mergeall.sh ${PN}*rescore.sc
# parse score files
makecsv.sh -merge -wt_label non-ts ${PN}.csv > ${PN}-in.csv
# run PSI-BLAST
run_psiblast.sh YBR109C
# make input file for classifier
prepare.sh ${PN}-in.csv
# make predictions
predict.sh ${PN}.arff
# view ranked predictions
less $PN-svmlin.txt $PN-svmrbf.txt
