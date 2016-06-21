for x in `ls -d */`; do 
    if [ -d $x ]; then
        cd $x
        x=`basename $x`
        if [ ! -f frags.${x}.25.pdb ]; then
            echo $x
        fi
        cd ../
    fi
done
