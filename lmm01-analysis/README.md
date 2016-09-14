This directory contains analysis of metagenomes mapped to complete *Microcystis* phage genomes

# Repo Organization: 

### Data:
Bam files for lmm01 and MamV from four lake erie viral metagenomes 

### R:
.Rmd and .html files with plot generation

### plots:
coverage and recruitment plots

### scripts: 
python, perl, bash scripts to get from bam files to R

# Data:

### Microcystis phage genomes
Microcystis phage genomes, lmm01 and MaMV, were pulled from genbank.    

### Samples
Running bowtie2 requires having the fastq files from each sample. These fastq files are stored on Greg's servers, and information on their location and sample numbers can be found [here](https://docs.google.com/spreadsheets/d/1MB-bntCNatELKqN2KJJiZtJ41XpM4rnIJbS4dBqwZ6Q/edit#gid=0).      

Directories were extracted from their tar archives with lines like this:     
``` nohup tar -xzvf /omics/HABs/Run_1351.tar.gz Run_1351/Sequences/Sample_49619/ Run_1351/Sequences/Sample_49618/ &```

Globus was used to transfer to flux.

# Protocol:
**1) Reads were mapped to each genome using bowtie2 with the following command:**       
```bowtie2 –no-unal –non-deterministic -p 20 -X 2000 -x path_to_db -1 R1_files.fastq.gz -2 R2_files.fastq.gz | samtools view -Sb - > Out_file.bam```

**Note: need to import script from flux**

Note: bowtie2 first tries to map the reads “consistently” (considering mate pair information), then separately if the pair doesn’t seem to match the same contig. This process is guided by the max insert size parameter (-X above). For our purposes, this was set at 2000. The insert size of the mapped Tara dataset centered around ~150 to 400bp, but anomolous much longer inserts did occur. Since this the `-X` is a maximum, an inflated (and unrealistic) value was used to cast a wide net and capture all mapped pairs. This has little to no effect on the mapping on itself, and mostly impacts the summary that bowtie2 provides (% of reads mapped consistently, inconsistently, etc).

The resulting bam files can be found in the directory **data** followed by  the subdirectory for each complete *Microcystis* genome. 

**2) Bam files need to be indexed and sorted.**       
The script `sort-index-bam.sh` is used for this purpose, but currently needs some updating to run generically.      
     
**3) bam_file_to_coverage.py transforms bam files to coverage.**       
This results in comma separated tables of coverage information. Each sample has a corresponding csv file that reports the number of bp mapped, number of reads mapped, corresponding coverage (the number of bp divided by genome length), and the % of the genome that is covered in each sample (i.e. `Sample_49639-lmm01-align.bam_sorted_coverage.csv`).    

After these are generated, a concatenated version needs to be made with all samples in one csv. Since working with a small number of samples, I did this in excel, but i should probably write a script. . . 

**4) mapping_to_ggplot_final.pl transforms bam files into a csv matrix to generate recruitment plots.**      
The resulting file is called `All_[genome]_recruitment.csv`. 

**5) Coverage and recruitment files are imported into R to generate plots.**    
Using the coverage csv file for all the samples, we plot the global coverage and % of genome mapped for each genome. Using the `All_[genome]_recruitment.csv` we generate a plot of the distribution of % of Id. Finally, recruitment plots for all viruses are generated, the results of which are output.

A composite plot is called `[genome]-combined-plot.pdf`.

Note: `Make_plots.R` is Simon's script which I used as the skeleton for the .Rmd file. Some important changes have been made, such as plotting the samples separately. 
