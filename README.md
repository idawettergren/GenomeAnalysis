# Genome Analysis

## Aim of project
The aim of this project is to reproduce some of the results from "The draft genome of tropical fruit durian (Durio zibethinus)" by Bin Tean Teh et al.

## Methods
A variety of analyses will be used during the project. These are Canu, Pilon and RepeatMaker for the DNA assembly, Trinity for the RNA assembly, BRAKER for annotation and BWA as well as STAR will be used as aligners.

## Timeline (from the Student Manual checkpoints)
In order to be on track the genome assembly should be done on April 15th, the annotation on April 29th and the RNA mapping on May 4th.

## Data organization
I don’t know how much space I will need to store my files, so the plan is to keep files compressed for as long as possible to conserve space. All files used in the project will be stored in my UPPMAX folder. This repository will contain all files besides the data. The files used in the project will be in many different file formats, for example fasta and fastq files.

# Log
A summary of what I have done.

## Retrieving metadata and raw data
Searched https://www.ncbi.nlm.nih.gov/sra for the accession number of the samples (PRJNA400310), with the accession number taken from the paper. Downloaded the runinfo table file ‘SraRunTable.txt’ that contain info, but no actual data, on all runs and move it to my metadata directory in UPPMAX.
Created soft links to raw pacbio and illumina data for scaffold 10 in the project folder from my raw_data folder (soft links to avoid using up a lot of space). The commands used for this can be found in 'retrieve_data'. 

## Canu genome assembly
Run canu genome assembly on the compressed pacbio raw data (‘SRR6037732_scaffold_10.fq.gz’) by running the script called ‘genome_assembly.sh’. First try submitted 03.28 April 13 because all my best work is done when I'm sleep deprived.




