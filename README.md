# Genome Analysis

## Aim of project
The aim of this project is to reproduce some of the results from "The draft genome of tropical fruit durian (Durio zibethinus)" by Bin Tean Teh et al. I will not analyze the entire genome, but will instead focus on only scaffold 10.

## Methods
A variety of analyses will be used during the project. These are Canu, Pilon and RepeatMaker for the DNA assembly, Trinity for the RNA assembly, BRAKER for annotation and BWA as well as STAR will be used as aligners.

## Timeline (from the Student Manual checkpoints)
In order to be on track the genome assembly should be done on April 15th, the annotation on April 29th and the RNA mapping on May 4th.

## Data organization
I don’t know how much space I will need to store my files, so the plan is to keep files compressed for as long as possible to conserve space. All files used in the project will be stored in my UPPMAX folder. This git repository will contain all files besides the data. The files used in the project will be in many different file formats, for example fasta and fastq files.
The files 'data_organisation.txt' and 'data_organisation.png' give an overview of how I have structured my files in UPPMAX.

# Log
An overview of what I have done so far.

## Retrieving metadata and raw data
I searched [NCBI](https://www.ncbi.nlm.nih.gov/sra) for the accession number of the samples (PRJNA400310), with the accession number taken from the paper. Downloaded the runinfo table file ‘SraRunTable.txt’ containing info on all runs and moved it to my metadata directory in UPPMAX.
Created soft links (to avoid using up a lot of space) to files with PacBio and Illumina reads for scaffold 10 in the project folder from my raw_datafolder. The commands used for this can be found in ['retrieve_data'](code/retrieve_data). 

## Canu genome assembly
Ran Canu genome assembly on the compressed PacBio reads (‘SRR6037732_scaffold_10.fq.gz’) by running the script [‘genome_assembly.sh’](code/genome_assembly.sh). In 'genome_assembly.sh' genomeSize is set to 24.2 because the data from [NCBI](https://www.ncbi.nlm.nih.gov/Traces/wgs/NSDW01?display=contigs), where scaffold 10's length is 24 162 007.
The first job was submitted 03.28 on April 13th, because working normal hours are for nerds.
First Canu assembly returned an error, the '-num_threads' parameter was invalid so I'm gonna use 'maxThreads' instead.
Second Canu job submitted 06.27 on April 13th and is currently running (hopefully without any problems).




