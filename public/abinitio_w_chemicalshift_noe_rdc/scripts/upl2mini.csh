cat $1 | grep -v ^\# | awk '($4-$1 > 3) || ($1-$4 > 3) {print}' | awk -v offset=$2 '(($4-offset)>0) && (($1-offset)>0) {printf("AmbiguousNMRDistance %5s %5d %5s %5d BOUNDED 1.5 %5.3f 0.3 NOE; rawdata %5.3f\n",$3,$1-offset,$6,$4-offset,$7+0.3/2,$7);}' 


