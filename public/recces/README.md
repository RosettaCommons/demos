# RECCES (Reweighting of Energy-function Collection with Conformational Ensemble Sampling)

KEYWORDS: RNA

## Author
This README was written in May 2015, by Fang-Chieh Chou (fcchou@stanford.edu).

## Brief
This demo illustrates the RECCES pipeline for computing the folding free energy of RNA helices. To keep it simple, we just show how to run RECCES simulated tempering (ST) on one construct.

Python codes needed to run the job are located at tools/recces/. You need to include this path into your PYTHONPATH to run the following demo. The python codes have only been tested with Python v2.7.

## Detail
We first do a quick a prerun. First create a separate folder:
```
    mkdir prerun
    cd prerun/
```
Now run the following commands:
```
    recces_turner -score:weights stepwise/rna/turner -rna:farna:thermal_sampling:n_cycle 300000 -seq1 gu -seq2 ac -rna:farna:thermal_sampling:temps 0.8 -rna:farna:thermal_sampling:out_prefix prerun
    recces_turner -score:weights stepwise/rna/turner -rna:farna:thermal_sampling:n_cycle 300000 -seq1 gu -seq2 ac -rna:farna:thermal_sampling:temps 1.0 -rna:farna:thermal_sampling:out_prefix prerun
    recces_turner -score:weights stepwise/rna/turner -rna:farna:thermal_sampling:n_cycle 300000 -seq1 gu -seq2 ac -rna:farna:thermal_sampling:temps 1.4 -rna:farna:thermal_sampling:out_prefix prerun
    recces_turner -score:weights stepwise/rna/turner -rna:farna:thermal_sampling:n_cycle 300000 -seq1 gu -seq2 ac -rna:farna:thermal_sampling:temps 1.8 -rna:farna:thermal_sampling:out_prefix prerun
    recces_turner -score:weights stepwise/rna/turner -rna:farna:thermal_sampling:n_cycle 300000 -seq1 gu -seq2 ac -rna:farna:thermal_sampling:temps 3.0 -rna:farna:thermal_sampling:out_prefix prerun
    recces_turner -score:weights stepwise/rna/turner -rna:farna:thermal_sampling:n_cycle 300000 -seq1 gu -seq2 ac -rna:farna:thermal_sampling:temps 7.0 -rna:farna:thermal_sampling:out_prefix prerun
    recces_turner -score:weights stepwise/rna/turner -rna:farna:thermal_sampling:n_cycle 300000 -seq1 gu -seq2 ac -rna:farna:thermal_sampling:temps 30  -rna:farna:thermal_sampling:out_prefix prerun
```
These preruns generate data for computing the ST weights in the following run. They may each take about 1 hour.

The ST weights can then be determined using the following code snippet:
```
    from recces.util import weight_evaluate
    weight_evaluate('./', 'prerun_hist_scores.gz')
```
Here it outputs a list of temperatures and the corresponding ST weights.

Now we create a new folder `ST` and run simulated tempering:
```
    mkdir ST
    recces_turner -score:weights stepwise/rna/turner -seq1 gu -seq2 ac -rna:farna:thermal_sampling:n_cycle 9000000 -rna:farna:thermal_sampling:temps 0.8 1 1.4 1.8 3 7 30 -st_weights 0 7.33 14.6 17.32 18.87 18.34 17.09 -rna:farna:thermal_sampling:out_prefix ST -save_score_terms
```
We also need to run a simulation at infinite temperature:
```
    recces_turner -score:weights stepwise/rna/turner -seq1 gu -seq2 ac -rna:farna:thermal_sampling:n_cycle 300000 -rna:farna:thermal_sampling:temps -1 -rna:farna:thermal_sampling:out_prefix kT_inf -save_score_terms
```
Note that here we use the "-save_score_terms" option to cache the contributions of each score term, so we may easily reweight the score function.

After the end of the run, the free energy can be computed using python codes:
```
    from recces.data import SingleSimulation, KT_IN_KCAL, N_SCORE_TERMS
    curr_wt = [0.73, 0.1, 0.0071, 0, 4.26, 2.46, 0.25, 0, 1.54, 4.54]
    sim = SingleSimulation('ST/', curr_wt)
    print sim.value, sim.value * KT_IN_KCAL
```
Here sim.value gives the free energy of the molecule in the unit of kT (T = 37C). To convert it to kcal/mol, we multiply the number by KT_IN_KCAL. In the `example_output` directory, the values in kT and in kcal/mol were: `1.35766728357 0.836772916126`. In an independent rerun, we got: `1.25722532755 0.774867389305`.

We may also reweight the score function and obtain the new free energy:
```
    import numpy as np
    sim.reweight(np.ones(N_SCORE_TERMS))
    print sim.value, sim.value * KT_IN_KCAL
```
In the `example_output` directory, we get: `31.8988957998 19.660289662`. In an independent rerun, we got: `31.9550773046 19.694916085`.

To determine the nearest neighbor energies, you need to add/subtract several free energies based on simulations of several single-stranded and double-stranded constructs.
