#!/bin/bash
protein=$1
outdir=pssm2
psiblast -query $protein.fasta -db nr -evalue 1 -inclusion_ethresh .001 -comp_based_stats 1 -outfmt 7 -num_descriptions 3000 -num_alignments 300 -num_threads 1 -num_iterations 2 -pseudocount 2 -out $outdir/$protein.pb -out_ascii_pssm $outdir/$protein.pssm -out_pssm $outdir/$protein.cp -export_search_strategy $outdir/$protein.ss
gzip -f $outdir/${protein}.{cp,pb}
