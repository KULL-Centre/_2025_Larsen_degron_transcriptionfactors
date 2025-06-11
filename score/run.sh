#!/bin/bash

# Merge counts files and map to libraries defined in ../library/transfac90_lib*.csv and control_dna.csv
#     get counts files from ERDA
Rscript merge_and_map.r samples.csv ../counts/*_counts.txt > merge_and_map.out

# Process FACS data
#     get FACS files from ERDA
cd facs
Rscript structure_facs.r data/*.csv > structure_facs.out
cd ..

# Calculate degron scores
Rscript scores.r counts.rda > scores.out
