#!/usr/bin/tcsh
if( $#argv < 1 ) then
  echo "enter pdb-id as arg-1" && exit -1
endif

#be sure to set these paths for the demo to work
set pathToDemo=/vol/ek/londonir/protocol_capture2/FlexPepDockAbInitio
set pathToVall=/vol/ek/share/rosetta/rosetta_fragments_latest/nnmake_database

set id=$1
set pepChain="B" #assumes peptide chain is B
ln -s $pathToDemo/input_files/2b1z.native.pdb ./native.pdb
ln -s $pathToDemo/input_files/2b1z.start.pdb  ./start.pdb
set nresA=`$pathToDemo/scripts/getFastaFromCoords.pl -pdbfile start.pdb | awk -F- '{print $2}'`
set nresB=`$pathToDemo/scripts/getFastaFromCoords.pl -pdbfile start.pdb -c B | awk -F- '{print $2}'`
cp $pathToDemo/input_files/flags .
cp $pathToDemo/input_files/prepack_flags .

# if peptide is shorter than 9 don't use 9-mer fragments
if ($nresB < 9) then
 sed -i 's/-flexPepDocking:frag9/#-flexPepDocking:frag9/g' flags
 sed -i 's/-frag9/#-frag9/g' flags
endif

# make fragments
mkdir frags 
cd frags
$pathToDemo/scripts/getFastaFromCoords.pl -pdbfile ../start.pdb -chain $pepChain > ${id}${pepChain}.fasta
ln -s ../native.pdb ref.pdb
ln -s $pathToVall/vall.dat.2006-05-05
cp $pathToDemo/input_files/frags/psi_L1.cfg .
sed "s/XXXX/$id/g" $pathToDemo/input_files/frags/flags > flags
if($nresB < 9) sed -i "s/9 5 3/5 3/g" flags
$pathToDemo/scripts/frags/make_fragments.pl ${id}${pepChain}.fasta -nosam > log.makeFrags
foreach i (3 5 9)
 if($i <= $nresB) then
   $pathToDemo/scripts/frags/shift.sh frags.${i}mers $nresA > frags.${i}mers.offset
 endif
end
cd ..
