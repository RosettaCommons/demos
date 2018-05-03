#!/bin/sh

#produces one pair of files, commandline and jobname. Commandline is used to execute jobname through lsf on the fleishman queue.

rm command & rm job/* &
rm H1_H2_like_length L1_L2_like_length L3_like_length H3_like_length segment_lengths_script_vars


tmp=$1
L1_L2_len=$(((tmp-1)*4))
tmp=$2
L3_len=$(((tmp+7)*4))
tmp=$3
H1_H2_len=$(((tmp-1)*4))
tmp=$4
H3_len=$(((tmp+9)*4))


while read p; do W=`echo $p |wc -w`; if [ "$W" == "$H1_H2_len" ]; then echo $p |awk '{print $NF}' >> H1_H2_like_length; fi; done < AB_db_files/H1_H2.db
while read p; do W=`echo $p |wc -w`; if [ "$W" == "$L1_L2_len" ]; then echo $p |awk '{print $NF}' >> L1_L2_like_length; fi; done < AB_db_files/L1_L2.db
while read p; do W=`echo $p |wc -w`; if [ "$W" == "$L3_len" ]; then echo $p |awk '{print $NF}' >> L3_like_length; fi; done < AB_db_files/L3.db
while read p; do W=`echo $p |wc -w`; if [ "$W" == "$H3_len" ]; then echo $p |awk '{print $NF}' >> H3_like_length; fi; done < AB_db_files/H3.db
	
for i in `seq 1 500`; do
	H3=`cat H3_like_length |sort --random-sort|head -1`;	
	L3=`cat L3_like_length |sort --random-sort|head -1`;
	H1_H2=`cat H1_H2_like_length|sort --random-sort|head -1`;
	L1_L2=`cat L1_L2_like_length|sort --random-sort|head -1`;
	


	name="$L1_L2"_"$L3"_"$H1_H2"_"$H3";
	#echo $name;
	echo -parser:script_vars entry_H1_H2="$H1_H2" entry_L1_L2="$L1_L2" entry_H3="$H3" entry_L3="$L3" >> segment_lengths_script_vars
	#./create_job.sh $name  @flags -parser:script_vars prefix=$name entry_H3=$H3 entry_H1_H2=$H1_H2 entry_L1_L2=$L1_L2 entry_L3=$L3 -out:prefix "$name"_ 
done
	

