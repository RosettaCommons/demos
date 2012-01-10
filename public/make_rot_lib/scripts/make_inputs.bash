#!/bin/bash

# set input template from command line
input_in=$1

# set phi ranges
phi_start=-180
phi_stop=180
phi_step=10
phi_string=XXX

# set psi ranges
psi_start=-180
psi_stop=180
psi_step=10
psi_string=YYY

# for each phi/psi pair
phi=$phi_start
while [ $phi -le $phi_stop ]; do
    psi=$psi_start
    while [ $psi -le $psi_stop ]; do

	    # make output filename
	    output_in=`echo $input_in | sed "s/$phi_string/$phi/g" | sed "s/$psi_string/$psi/g"`

	    # copy inputs to outputs
	    cp $input_in $output_in

	    # do replacements on output
	    sed -i "" "s/$phi_string/$phi/g" $output_in
	    sed -i "" "s/$psi_string/$psi/g" $output_in

	    # output to let user know what we are up to 
	    echo $phi $psi $output_in

	psi=$(( $psi + $psi_step ))
    done
    phi=$(( $phi + $phi_step ))
done
