#!/usr/bin/tcsh
set PATH_TO_EXE=/vol/ek/ravehb/rosetta/svn_mini/rosetta_source/bin/
set PATH_TO_DB=/vol/ek/ravehb/rosetta/svn_mini/rosetta_database/
set PATH_TO_DEMO=/vol/ek/ravehb/rosetta/RosettaCon2011/demos/global_dock_ssrA_peptide_against_sspB
set PATH_TO_SCRIPTS=$PATH_TO_DEMO/scripts/
set PATH_TO_UTILS=$PATH_TO_SCRIPTS/utils/
set PATH_TO_VALL=/vol/ek/share/rosetta/rosetta_fragments_latest/nnmake_database
set pdbid=1OU8
set origRecepChain="B"
set origPepChain="D"

# ==========================
# create folder
cd input_files
cp ../starting_files/$pdbid.pdb.gz .
gunzip $pdbid.pdb.gz

# prepare clean native complex file, with receptor as chain A and peptide as chain B
if(! -e native.pdb) then
 $PATH_TO_UTILS/extract_chains_and_range.pl -p ${pdbid}.pdb -c $origRecepChain -o ${pdbid}_native_receptor.pdb
 $PATH_TO_UTILS/switchChain.pl ${pdbid}_native_receptor.pdb $origRecepChain A
 $PATH_TO_UTILS/extract_chains_and_range.pl -p $pdbid.pdb -c $origPepChain -o ${pdbid}_native_peptide.pdb
 $PATH_TO_UTILS/switchChain.pl ${pdbid}_native_peptide.pdb $origPepChain B
 cat ${pdbid}_native_receptor.pdb > ${pdbid}_native_complex.pdb
 echo TER >> ${pdbid}_native_complex.pdb
 cat ${pdbid}_native_peptide.pdb >> ${pdbid}_native_complex.pdb
 echo TER >> ${pdbid}_native_complex.pdb
 ln -s ${pdbid}_native_complex.pdb native.pdb
 echo native.pdb created
endif

#create complex between native receptor and extended peptide from FASTA sequence
if(! -e start.pdb) then
 $PATH_TO_UTILS/getFastaFromCoords.pl -pdbfile ${pdbid}_native_peptide.pdb -c B > peptide.fasta
 $PATH_TO_EXE/BuildPeptide.linuxgccrelease -in:file:fasta peptide.fasta -database $PATH_TO_DB -out:file:o ${pdbid}_peptide_extended_from_fasta.pdb
 $PATH_TO_UTILS/switchChain.pl ${pdbid}_peptide_extended_from_fasta.pdb A B
 cat ${pdbid}_native_receptor.pdb > ${pdbid}_complex_with_extended_peptide.pdb
 echo TER >> ${pdbid}_complex_with_extended_peptide.pdb
 cat ${pdbid}_peptide_extended_from_fasta.pdb >> ${pdbid}_complex_with_extended_peptide.pdb
 echo TER >> ${pdbid}_complex_with_extended_peptide.pdb
 ln -s ${pdbid}_complex_with_extended_peptide.pdb start.pdb
 echo start.pdb created
endif

#prepack complex with extended peptide
if(! -e start.ppk.pdb) then
 $PATH_TO_EXE/FlexPepDocking.linuxgccrelease -s start.pdb -native native.pdb \
 -database $PATH_TO_DB -scorefile score.ppk.sc @prepack_flags > log.prepack
 mv start_0001.pdb start.ppk.pdb
 echo start.ppk.pdb prepacked model created
endif

# compute nres
set nresA=`awk '/^ATOM/ && substr($0,14,3)=="CA " && substr($0,22,1)=="A"' start.pdb | grep "CA" | wc -l`
set nresB=`awk '/^ATOM/ && substr($0,14,3)=="CA " && substr($0,22,1)=="B"' start.pdb | grep "CA" | wc -l`
echo $nresA residues detected in chain A
echo $nresB residues detected in chain B 

# make fragments
cd frags
$PATH_TO_UTILS/getFastaFromCoords.pl -pdbfile ../start.pdb -chain B > ${pdbid}B.fasta
$PATH_TO_UTILS/extract_chains_and_range.pl -p ../native.pdb -c B -o ref.pdb
ln -s $PATH_TO_VALL/vall.dat.2006-05-05
sed -i "s/XXXX/$pdbid/g" flags 
if($nresB < 9) sed -i "s/9 5 3/5 3/g" flags
$PATH_TO_SCRIPTS/frags/make_fragments.pl ${pdbid}B.fasta -nosam > log.makeFrags
foreach i (3 5 9)
 if($i <= $nresB) then
   $PATH_TO_SCRIPTS/frags/shift.sh frags.${i}mers $nresA > frags.${i}mers.offset
   if(-s frags.${i}mers.offset) then
    echo "Fragments file created!"
   else
    echo "ERROR: file frags.${i}mers.offset is either empty or not created"
   endif
 endif
end
cd ..


