#!/usr/bin/awk -f

BEGIN {
RS="\n"

# initialize the atom counter and residue counter
natom = 0
nresidue = 0
cur_res=-10

}

/REMARK/ {
print $0 > "renum.pdb"
next
}

/HETATM/ {
next
}

/CONECT/ {
next
}

/ATOM/ {

if($1 != "ATOM") next

# get our new atom and residue numbers
natom++

if(cur_res != $6)
        {
        cur_res=$6
        nresidue++
        }

# paste in the atom number

frnt_str = substr($0,0,6)
back_str = substr($0,12,75)

if(natom < 10)
        workstring = sprintf("%s    %d%s",frnt_str,natom,back_str)
else if(natom < 100)
        workstring = sprintf("%s   %d%s",frnt_str,natom,back_str)
else if(natom < 1000)
        workstring = sprintf("%s  %d%s",frnt_str,natom,back_str)
else if(natom < 10000)
        workstring = sprintf("%s %d%s",frnt_str,natom,back_str)
else
        workstring = sprintf("%s%d%s",frnt_str,natom,back_str)

#paste in the residue number

frnt_str = substr(workstring,0,21)
back_str = substr(workstring,28,65)

if(nresidue < 10)
        workstring = sprintf("%s    %d %s",frnt_str,nresidue,back_str)
else if(nresidue < 100)
        workstring = sprintf("%s   %d %s",frnt_str,nresidue,back_str)
else if(nresidue < 1000)
        workstring = sprintf("%s  %d %s",frnt_str,nresidue,back_str)
else if(nresidue < 10000)
        workstring = sprintf("%s %d %s",frnt_str,nresidue,back_str)
else
        workstring = sprintf("%s%d %s",frnt_str,nresidue,back_str)

# get the first 72 characters
tempstring=substr(workstring,0,76)

# append the appropriate string for our segment
newtempstring = sprintf("%s%s",tempstring,strrange[counter])

print newtempstring > "renum.pdb"

next

}

{
print > "renum.pdb"
}

