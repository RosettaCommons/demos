###############################################################################################
###                                                                                         ###
###       Rosetta command lines for the design protocol                                     ###
###       ( Compatible with the Rosetta SVN code repository after the Revision 51602 )      ###
###                                                                                         ###
###       Principles for designing ideal protein structures                                 ###
###       Nobuyasu Koga, Rie Tatsumi-Koga, Gaohua Liu, Rong Xiao, Thomas B. Acton,          ###
###       Gaetano T. Montelione, and David Baker                                            ###
###       Nature 2012, Supplementary Data                                                   ###
###                                                                                         ###
###############################################################################################

# Run the following command line
~/rosetta/rosetta_source/bin/rosetta_scripts.linuxgccrelease -database ~/rosetta/rosetta_database -s ./input.pdb -parser:protocol ./input.xml -nstruct 100 -ex1 -ex2aro
