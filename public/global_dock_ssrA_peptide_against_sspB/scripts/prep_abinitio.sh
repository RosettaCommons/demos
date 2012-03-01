#!/usr/bin/tcsh
#be sure to set these correctly
set pdbid=1ou8
set recepChain="A"
set pepChain="B" #assumes peptide chain is B
set pathToVall=/vol/ek/share/rosetta/rosetta_fragments_latest/nnmake_database
set pathToDemo=/vol/ek/ravehb/rosetta/RosettaCon2011/demos/global_dock_ssrA_peptide_against_sspB/
set pathToUtils=$pathToDemo/scripts/utils/
set PATH_TO_ROSETTABIN=/vol/ek/ravehb/rosetta/svn_mini/rosetta_source/bin/
set PATH_TO_ROSETTADB=/vol/ek/ravehb/rosetta/svn_mini/rosetta_database/

# =======================
# assumes model is already in start_model
cd $pathToDemo/output
ln -s ./start_model/start.ppk.pdb
ln -s ./start_model/native.pdb
cp $pathToDemo/input_files/flags .
ln -s $pathToDemo/input_files/site_constraints.cst # constraints file specifying peptide binding site for low-res ab-initio step

set nresA=`awk '/^ATOM/ && substr($0,14,3)=="CA " && substr($0,22,1)=="'"$recepChain"'"' start.pdb | grep "CA" | wc -l`
set nresB=`awk '/^ATOM/ && substr($0,14,3)=="CA " && substr($0,22,1)=="'"$pepChain"'"' start.pdb | grep "CA" | wc -l`
echo $nresA residues detected in chain $recepChain 
echo $nresB residues detected in chain $pepChain 

# if peptide is shorter than 9 don't use 9-mer fragments
if ($nresB < 9) then
 sed -i 's/-flexPepDocking:frag9/#-flexPepDocking:frag9/g' flags
 sed -i 's/-frag9/#-frag9/g' flags
endif

# make fragments
mkdir frags 
cd frags
$pathToUtils/getFastaFromCoords.pl -pdbfile ../start.pdb -chain $pepChain > ${pdbid}${pepChain}.fasta
$pathToUtils/extract_chains_and_range -p ../native.pdb -c $pepChain -o ref.pdb
ln -s $pathToVall/vall.dat.2006-05-05
cp $pathToDemo/input_files/frags/psi_L1.cfg .
sed "s/XXXX/$pdbid/g" $pathToDemo/input_files/frags/flags > flags
if($nresB < 9) sed -i "s/9 5 3/5 3/g" flags
$pathToDemo/scripts/frags/make_fragments.pl ${pdbid}${pepChain}.fasta -nosam > log.makeFrags
foreach i (3 5 9)
 if($i <= $nresB) then
   $pathToDemo/scripts/frags/shift.sh frags.${i}mers $nresA > frags.${i}mers.offset
 endif
end
cd ..
echo "Fragments file created!"
echo

# prepack structure
$PATH_TO_ROSETTABIN/FlexPepDocking.linuxgccrelease -database $PATH_TO_ROSETTADB \
    -s start.pdb -native native.pdb @prepack_flags > prepack.log
mv start_0001.pdb start.ppk.pdb

echo "Pre-packing done!"
echo
