#!/bin/bash
#
# 
# input file structure
# pdb_file    chainID    start end

func()
{
database="/path/to/rosetta/rosetta_database"
window=5
echo "Input file: $1"
wd=`pwd`
mod1=1
mod2=2
mod3=3
mod4=0
j=1
for i in `cat $1`; do
    let "k= $(($j % 4))" #allow getting every 4th input
    if [ $k -eq $mod1 ]
	then
	file=$i
    fi
    if [ $k -eq $mod2 ]
	then
	chain=$i
    fi
    if [ $k -eq $mod3 ]
	then
	start=$i
    fi
    if [ $k -eq $mod4 ]
	then
	end=$i
#runs the command now
	cmd="./homodimer_maker.release -database $database -s $file \
-run::chain $chain -sheet_start $start -sheet_stop $end -window_size $window \
-ignore_unrecognized_res -mute core protocols.moves.RigidBodyMover "

	echo $cmd
	$cmd
    fi


#echo $i
let "j += 1"
done
}

# clean up old directories, make fresh ones
#rm tracer_output.log
func $1
echo "Done."