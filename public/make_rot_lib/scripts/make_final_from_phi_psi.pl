#!/usr/bin/perl -w
# Make a set of input files for the make_rot_lib protocol

@phi = ( -180, -170, -160, -150, -140, -130, -120, -110, -100, -90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180);
@psi = ( -180, -170, -160, -150, -140, -130, -120, -110, -100, -90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180);

@aa = ("UNK");

$final_filename = "$ARGV[0].rotlib";

foreach $aa (@aa) {

    open(FINAL, "> $final_filename") or die "Cant open file $final_filename!!!\n";

    foreach $phi ( @phi ) {
	foreach $psi ( @psi ) {
	    #create filename and open file
	    $filename = $aa . "_" . $phi . "_" . $psi . ".rotlib";
	    print $filename . "\n";
	    open(FILE, "$filename") or die "Cant open file $filename!!!\n";
	    @lines = <FILE>;
	    
	    #outptut data to file
	    print FINAL @lines;
	}
    }
}
