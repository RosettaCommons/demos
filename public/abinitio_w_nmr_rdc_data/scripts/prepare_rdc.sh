#input:
more $1 |awk '{if (NR>1) {print $1," H ",$1," N ",$NF}}'>med1.rdc
