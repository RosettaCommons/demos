#!/usr/bin/tcsh

sed 's/ position/position/g' $1 | awk '{ if ( substr ( $0,1,3 ) == "pos" ) {print substr ( $0,0,18 ) sprintf ("%4d",substr ( $0,19,4 ) + '"$2"' ) substr ( $0,23,1000 ) ; } else {print ; }}'
