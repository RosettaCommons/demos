#/bin/bash
set -e

rosetta_path=/nas02/home/k/h/khouli/rosetta
dbpath=${rosetta_path}/rosetta_database

basename=`pwd`
jobname=$1
struct="$basename""/""$2"
scriptfile="$basename""/""$3"
nstruct=100
trials=20
localname="$jobname.script.xml"

numprocs=128
que="day"
exepath="${rosetta_path}/rosetta_source/bin/rosetta_scripts.mpi.linuxgccrelease"

logfile="$jobname.run.log"

bsubcommand="bsub -q $que -n $numprocs -J $jobname -o $logfile.%J -a mvapich mpirun"

cmd="$bsubcommand $exepath -database $dbpath \
-l $struct \
-nstruct $nstruct \
-jd2:ntrials $trials \
-out:pdb_gz true \
-parser:protocol $localname \
-mpi_tracer_to_file $jobname.tracers \
-overwrite \
-chemical:exclude_patches LowerDNA UpperDNA Cterm_amidation SpecialRotamer protein_cutpoint_upper protein_cutpoint_lower VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm \
-score:weights score12minpack \
-atomic_burial_cutoff 0.01 \
-sasa_calculator_probe_radius 1.2 \
-ignore_unrecognized_res \
-no_optH false \
-skip_set_reasonable_fold_tree \
-no_his_his_pairE \
-fa_max_dis 6.0 \
-options:user \
"

[ -e "${jobname}/" ] && rm -r "${jobname}/"
mkdir $jobname
cd $jobname
cp $scriptfile $localname 
cp $struct start_struct.list
echo $cmd > command_run
$cmd
