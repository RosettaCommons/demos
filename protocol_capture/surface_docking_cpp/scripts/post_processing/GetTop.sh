#!/bin/bash

prefix=$1
num=$2

grep 'Total weighted' $prefix*.pdb | sort -nk5 | awk '{print $1, $5}' | sed s/://g > $prefix.sorted

mkdir TOP$num.$prefix

i=1
while [ $i -le $num ]
  do
  file=$( head -$i $prefix.sorted | tail -1 | awk '{print $1}' )
  echo Copying $file 
  cp $file TOP$num.$prefix
  i=$[ $i + 1 ]

done
