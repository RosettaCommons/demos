bcl svn revision: 4351

A plot can be made with angles of chi1 and chi2 on the x and y axes. The plot
shows the frequency with which a given chi1/chi2 angle pair is sampled by
rosetta for a given mtssl mutant. The dark black box shows the native chi1 and
chi2 angles +- 30 degrees. The light grey boxes show the chi1/chi2 angles that
are contained in the rotamer library. In order to make a plot for a set of
models of a single mutant of lysozyme, the command below can be used.
bin/bcl.exe AnalyzeRestraintAgreement -pdb_list input/pdbs.ls -analysis_prefix chi_angle_pair_distr_ -analysis_type_enumerated 'ChiAnglePairDistribution ( filename_postfix=.ChiAnglePairDistribution,histogram_binsize=30,collector_type=CollectorAAType ( aa_list ( methanesulfonothioate ) ) ,chi_angle_a=e_One, chi_angle_b=e_Two,set_title_and_label=1,title=A082,show_color_box=0,font_size=10,x_pixels=350,y_pixels=350,reference_chi_filename=input/chi_boxes.txt,font=/usr/share/fonts/dejavu-lgc/DejaVuLGCSansMono.ttf, grey_scale=1, plot_ratio=nan, bins_per_tic=1, center_tics=0,min_z=0.0,max_z=0.5,normalize_histogram=1) ' -message_level Standard > & chi_angle_pair_distribution.log &
To create the plot png 
gnuplot chi_angle_pair_distr_.ChiAnglePairDistribution
will create 
AnalyzeChiAnglePairDistribution.gnuplot.png

To get the percentage of models that agree with the native conformation, the
below command line can be used. It calculates the percentage and counts of how
frequently each proceeding chi angle is correctly recovered. For example, the
recovery of chi3 is the rate of correctly recovering chi1 and chi2 and chi3.
bin/bcl.exe AnalyzeRestraintAgreement -pdb_list input/pdbs.ls 0 -analysis_type_enumerated 'ChiAngleRecovery ( filename_postfix=.ChiAngleRecovery,collector_type=CollectorAAType ( aa_list ( methanesulfonothioate ) ) ,reference_chi_filename=input/native_chi.txt,tolerance=30,angle_unit=degree ) ' -analysis_prefix out_chi_recovery >& chi_angle_recovery.log & 
The generated text file contains the calculated numbers as well as the list of
protein models which maximally fulfill the native chi angles. If there are no
models which recover all 5 chi angles, the list of models will be the models
which recover as many chi angles as possible. The generated file is
out_chi_recovery.ChiAngleRecovery

