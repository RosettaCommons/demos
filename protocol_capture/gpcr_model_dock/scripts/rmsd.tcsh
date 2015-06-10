#!/bin/tcsh

foreach f (`ls *.pdb`)

  set name = `basename $f .pdb`

  grep HETATM $f > $name.het

  sed -i 's/HETATM/ATOM  /g' $name.het

  babel -ipdb $name.het -osd sdf.tmp

  cat sdf.tmp >> all.sdf

end

bzip2 all.sdf

rm *.tmp *.het


