bcl svn revision: 4434

Finding the subensemble of models that fits with the experimental distribution
can be accomplished using the command line provided below.
bin/bcl.exe FitEPRDistribution -exp_hist_single input/t4l_epr_059_A_159_A.histogram -exp_hist_data_columns 0 1 -model_list input/pdbs.ls 0 -start_size_range 5 20 -num_fits 1 -terminate_criteria 0.1 2500 -use_pdbid_numbering -message_level Standard -random_seed 12345 -convert_nanometers_to_angstrom -full_atom_sl_score -aaclass AAComplete -prefix fit_ >& fit.log &
This generates the files below which contain, respectively. 1) the distance distribution histogram of the resulting ensemble 2) the list of models that make up the resulting ensemble 3) the score of the ensemble distribution with the experimental distribution
fit_ensemble_00000.histogram  
fit_ensemble_00000.ls
fit_ensemble_00000.sc

To see a heatmap of the agreement of the fitted subensemble models' distribution versus the epr distribution the command line below can be used. The experimental mean and standard deviation is provided in the command line as a gaussian function.
bin/bcl.exe AnalyzeRestraintAgreement -pdb_list fit_ensemble_00000.ls -analysis_prefix coord_dstnc_distr_ -analysis_type_enumerated 'CoordinateDistanceDistribution(filename_postfix=.CoordinateDistanceDistribution,coord_a=LocatorCoordinatesAverage(LocatorAtom(locator_aa=LocatorAA(locator_chain=A,seq_id=59,use_pdb_id=1),atom_type=O1),LocatorAtom(locator_aa=LocatorAA(locator_chain=A,seq_id=59,use_pdb_id=1),atom_type=N1)),coord_b=LocatorCoordinatesAverage(LocatorAtom(locator_aa=LocatorAA(locator_chain=A,seq_id=159,use_pdb_id=1),atom_type=O1),LocatorAtom(locator_aa=LocatorAA(locator_chain=A,seq_id=159,use_pdb_id=1),atom_type=N1)),function=GaussianFunction(mean=41.9,sigma=2.7),histogram_minimum=0,histogram_binsize=2,histogram_num_bins=30,title="Mutant A 059 B 159",pixel_x=400,pixel_y=100,font="/usr/share/fonts/dejavu-lgc/DejaVuLGCSansMono.ttf",font_size=8,grey_scale=1,plot_ratio=0.25, mean_stddev_out_file=coord_dstnc_distr_ensemble.meanstdev.txt,bins_per_tic=2,center_tics=0,min_z=0,max_z=0.5,normalize_histogram=1)' -message_level Critical > & coord_dstnc_distr.log &
This generates two files. The first is coord_dstnc_distr_ensemble.meanstdev.txt which contains the mean and standard deviation of the distances between spin labels observed in the ensemble of models. The second file is a gnuplot script that can be used to generate the heat map. To generate the heatmap:
gnuplot coord_dstnc_distr_.CoordinateDistanceDistribution
This generates the file heatmap.png

To see a line plot of the agreement of the fitted subensemble models'
distribution versus the epr distribution the command line below can be used.
bin/bcl.exe FitEPRDistribution -distributions_to_lineplots fit_ensemble_00000.histogram line_plot_outfile.gnuplot true 0 80 -message_level Debug >& lineplot.log &
Then to generate the line plot:
gnuplot line_plot_outfile.gnuplot

A line plot of the integral of the probability distributions can provide a
complementary view of the agreement and can be obtained with the command line
below:
bin/bcl.exe FitEPRDistribution -distributions_to_lineplots_sum fit_ensemble_00000.histogram lineplot_sum.gnuplot 0 80 -message_level Debug > & lineplot_sum.log &
To generate the plot:
gnuplot lineplot_sum.gnuplot
