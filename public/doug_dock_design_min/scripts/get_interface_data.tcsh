#!/bin/tcsh

# get a list of pdb file names for join to work properly
gunzip *.gz
ls -l *.pdb | awk '{print $9 '} > out.TEMP1

echo "BEFORE LOOP"
head -n1 out.TEMP1

# get interface analyzier info from all the pdbs
 foreach i ( ENERGY SASA HB_ENER PACK )
     foreach j ( COMPLEX SEPERATE DIFF )
	echo ${i}_${j}
 	grep ${i}_${j} *.pdb > out.TEMP2
	sed -i.bak 's/\:/      /g' out.TEMP2
	echo "BEFORE JOIN"
	head -n1 out.TEMP1
 	join out.TEMP1 out.TEMP2 > out.TEMP3
	mv out.TEMP3 out.TEMP1
	echo "AFTER JOIN"
	head -n1 out.TEMP1
     end
end

cp out.TEMP1 out.TEMP_NO_SEQ

## get sequence info 
#echo $1
#grep $1 *.pdb  > out.TEMP2
#sed -i 's/\:/      /g' out.TEMP2
#awk '{ print $1 " " $2 '} out.TEMP2 > out.TEMP2A
#join out.TEMP1 out.TEMP2A > out.TEMP3
#mv out.TEMP3 out.TEMP1
#
#echo $2
#grep $2 *.pdb  > out.TEMP2
#sed -i 's/\:/      /g' out.TEMP2
#awk '{ print $1 " " $2 '} out.TEMP2 > out.TEMP2B
#join out.TEMP1 out.TEMP2B > out.TEMP3
#mv out.TEMP3 out.TEMP1
#
#echo $3
#grep $3 *.pdb  > out.TEMP2
#sed -i 's/\:/      /g' out.TEMP2
#awk '{ print $1 " " $2 '} out.TEMP2 > out.TEMP2C
#join out.TEMP1 out.TEMP2C > out.TEMP3
#mv out.TEMP3 out.TEMP1
#
#echo $4
#grep $4 *.pdb  > out.TEMP2
#sed -i 's/\:/      /g' out.TEMP2
#awk '{ print $1 " " $2 '} out.TEMP2 > out.TEMP2D
#join out.TEMP1 out.TEMP2D > out.TEMP3
#mv out.TEMP3 out.TEMP1
#
#echo $5
#grep $5 *.pdb  > out.TEMP2
#sed -i 's/\:/      /g' out.TEMP2
#awk '{ print $1 " " $2 '} out.TEMP2 > out.TEMP2E
#join out.TEMP1 out.TEMP2E > out.TEMP3
#mv out.TEMP3 out.TEMP1


# print it all out 
#head -n1 out.TEMP1 | awk '{ print "PDB\t" $2 "\t" $4 "\t" $6 "\t" $8 "\t" $10 "\t" $12 "\t" $14 "\t" $16 "\t" $18 "\t" $20 "\t" $22 "\t" $24 "\tSEQ1\tSEQ2\tSEQ3\tSEQ4\t" '} > out.TEMP4
#cat out.TEMP1 | awk '{ print $1 "\t" $3 "\t" $5 "\t" $7 "\t" $9 "\t" $11 "\t" $13 "\t" $15 "\t" $17 "\t" $19 "\t" $21 "\t" $23 "\t" $25 "\t" $26 "\t" $27 "\t" $28 "\t" $29 "\t" $30'} >> out.TEMP4
cat out.TEMP1 | awk '{ print $1 "\t" $3 "\t" $5 "\t" $7 "\t" $9 "\t" $11 "\t" $13 "\t" $15 "\t" $17 "\t" $19 "\t" $21 "\t" $23 "\t" $25'} >> out.TEMP4
sort -nuk 4 out.TEMP4 > out.TEMP5
mv out.TEMP5 out.ALL
rm out.TEMP*
