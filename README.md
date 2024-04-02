# Genes_from_histone_mod_122
The goal is to identify genomic intervals overlapping genes based on histone modification data.

This project uses processed data from a Chromatin Immunoprecipitation (ChIP-
seq) experiment for two different histone modifications from the same cell type. Both histone
modifications have a broad profile and are typically found within annotated gene bodies.
The processed data is for consecutive 200bp intervals of a 100 MB sequence portion of one
chromosome. The data has been processed such that each interval is assigned one of the
following four characters: {x,y,z,n}, which have the following meaning:
x – the number of ChIP-seq reads for only mark X in the interval was above a significance
threshold
y – the number of ChIP-seq reads for only mark Y in the interval was above a significance
threshold
z – the number of ChIP-seq reads for both marks X and Y in the interval was above a
significance threshold
n – the number of ChIP-seq reads for neither marks X and Y in the interval was above a
significance threshold
The task is to predict 200bp intervals that are most likely overlapping annotated gene bodies of
protein coding genes.
Note you are not given a labeled dataset for training. For the purposes of this project, the
identity of the specific ChIP-seq experiment and what portion of the genome the provided data
is from will be kept anonymous.

Inputs:
Given a single input file called input.fasta. Following a header line, the file contains
500,000 lines each containing one character in the set {n,x,y,z} as shown below. Each line
corresponds to a consecutive 200bp interval.
>seq
n
n
x
x
z
z
y
y
...
Project Output: The output for this project is your top 50,000 predictions of 200bp intervals that
overlap with the body of an annotated protein coding gene. You will provide the line numbers of
them one per line. Lines are numbered starting from 1. Below is a portion of a sample output
where among the top 50,000 predictions are lines corresponding to the 1st , 2nd, and 4th 200bp
interval.
1
2
4
...
