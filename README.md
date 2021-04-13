# Genome Analysis

## Aim of project
The aim of this project is to reproduce some of the results from ["The draft genome of tropical fruit durian (Durio zibethinus)" by Bin Tean Teh et al](https://www.nature.com/articles/ng.3972/). I will not analyse the entire genome, but will instead focus only on scaffold 10.

## Methods
A variety of analyses will be used during the project. These are Canu, Pilon and RepeatMaker for the DNA assembly, Trinity for the RNA assembly, BRAKER for annotation and BWA as well as STAR will be used as aligners.

## Timeline (from the Student Manual checkpoints)
In order to be on track the genome assembly should be done on April 15th, the annotation on April 29th and the RNA mapping on May 4th.

## Data organization
The files [data_organisation.txt](data_organisation.txt) and [data_organisation.png](data_organisation.png) give a simple overview of how I have organised my files.
I don’t know how much space I will need, so the plan is to keep files compressed for as long as possible and to use soft links for the raw data to conserve space. All files used in the project will be stored in my UPPMAX directory. This git repository will contain all files besides the data. The files used in the project will be in many different file formats, for example fasta and fastq files.

# Log
An overview of what I have done so far.

## Retrieving metadata and raw data
I searched [NCBI](https://www.ncbi.nlm.nih.gov/sra) for the accession number of the samples (PRJNA400310), with the accession number taken from the paper. Then I downloaded the runinfo table file ‘SraRunTable.txt’ containing info on all runs and moved it to my metadata directory.
Since the files with the raw data reads are big and I'm limited on storage space, I created soft links to the PacBio and Illumina reads for scaffold 10 in the project directory from my raw_data directory instead of copying the files. The commands used for creating the soft links can be found in [retrieve_data](code/retrieve_data). 

## Canu genome assembly
I ran a Canu genome assembly on the compressed PacBio reads in 'SRR6037732_scaffold_10.fq.gz' by running the script [genome_assembly.sh](code/genome_assembly.sh). In the script the parameter 'genomeSize' is set to 24.2m because of data from [NCBI](https://www.ncbi.nlm.nih.gov/Traces/wgs/NSDW01?display=contigs), where scaffold 10's length is stated to be 24 162 007.

The first Canu assembly job was submitted 03.28 on April 13th and returned an error as soon as the job started. The [output file](analysis/genome_assembly/canu/01_canu_assembly/slurm-4411986.out) states that the -num_threads parameter was invalid, so I will use the maxThreads parameter instead for the next job. A parameter that specifies the number of threads the job can use is needed in order for the job to use all the requested cores since without it the job will still be able to run on one core, which is a waste of data as well as the job probably not being able to finish in time.

The second Canu job was submitted 06.27 on April 13th and *didn't* return an error right away.
I was however nervous about either accidentaly requesting too many cores or not using all the requested cores, so I generated a graph over CPU usage for the job using the jobinfo command (full command with parameters can be found in [misc](code/misc)). The first graph of the job was made after 1h, 43min and 23s and can be found in [job_info_canu_01.png](job_info/job_info_canu_01.png) together with a short explanatory [text file](job_info/job_info_canu_01.txt). The job hadn't been running for too long, but the blue line (core usage) looks to be consistently close to the max capacity of 400% core usage, which would mean that the job is fully using all 4 requested cores.

# To-do list
* Update project plan according to feedback
* Pre-processing and trimming of Illumina reads
* Run the second Canu job (should be done around 23.30 on April 13th if the estimated 17h runtime is correct)
* Make the data_organisation image not look like shit
* Add, commit and push output files from the second Canu job to git
