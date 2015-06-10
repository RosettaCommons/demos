set terminal postscript eps enhanced color "Helvetica" 16
set output "example_outputsNGK_1a8d_score_vs_RMSD.eps"
set title "example-outputsNGK-1a8d final centroid and fa models"
set xlabel "local backbone RMSD"
set xtics out
set xrange [0:*]
set ylabel "Rosetta Energy Units"
set ytics out
set ytics nomirror
plot  "example_outputsNGK_1a8d_fa_score_vs_RMSD.dat" using 2:1 pt 6 lc rgb "red" title "full-atom" axes x1y1, \
 "example_outputsNGK_1a8d_fa_score_vs_RMSD.dat" index 1:1 using 2:1 pt 13 ps 2 lc rgb "orange" title "fa top10" axes x1y1, \
 "example_outputsNGK_1a8d_fa_score_vs_RMSD.dat" index 2:2 using 2:1 pt 12 ps 2 lc rgb "black" title "fa top-scoring" axes x1y1;

