#!/bin/sh

#grep out lines containing "total energy" in the pdbs, then sort by the second column, which will be the energies.
grep -H total_energy 1FNA_del__*pdb | sort -n -k2 