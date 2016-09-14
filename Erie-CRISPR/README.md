# data
The data section is organized with directories for the 14 bacterial metagenomes. To reduce running time for CRASS, we extracted two microcystis DR
types from each bacterial metagenome using grep. We then ran CRASS on those grepped reads. Within each sample directory, there are fasta
files with the spacers for the DR types. There are also additional directories containing the output of running blast between the spacers and 
viral metagenomes from samples 0024, 0048, 0108, 0132. 

# analysis
Rmd, md, and html files for R analysis done on CRASS output. The first part of the analysis focuses on basic summary statistics of CRISPRS.
The second part examines blast hits between spacers and the viral metagenomes. 

# Figs
These are figures automatically extracted from the R script in analysis
