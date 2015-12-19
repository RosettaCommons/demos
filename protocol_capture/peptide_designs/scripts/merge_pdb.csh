# csh to merge two pdb files

foreach file (`cat $argv[2]`)
    head -5432 $argv[1] > $file.m.pdb
    tail -426 $file.pdb >> $file.m.pdb
end
