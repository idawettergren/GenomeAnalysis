# Genome Analysis

## Aim of project
The aim of this project is to reproduce some of the results from ["The draft genome of tropical fruit durian (Durio zibethinus)" by Bin Tean Teh et al](https://www.nature.com/articles/ng.3972/). I will not analyse the entire genome, but will instead focus only on scaffold 10.

## Methods (incomplete)
A variety of analyses will be used during the project. To assemble the PacBio reads I will use Canu, for quality control of the Illumina reads I will use FastQC and thereafter trim them using Trimmomatic. Then the assembled PacBio reads can be corrected with the help of the trimmed Illumina reads using BWA, after which the assembly can be further improved with Pilon and RepeatMasker. The assembly can then be annotated using BRAKER.

## Timeline (from the Student Manual checkpoints)
In order to be on track the genome assembly should be done on April 15th, the annotation on April 29th and the RNA mapping on May 4th.

## Data organization
The files [data_organisation.txt](data_organisation.txt) and [data_organisation.png](data_organisation.png) give a simple overview of how I have organised my files.
I don’t know how much space I will need so the plan is to keep files compressed for as long as possible and to use soft links for raw data to conserve space. All files used in the project will be stored in my UPPMAX directory.
This git repository will contain all files except for the data.
The files used in the project will be in many different file formats, for example fasta and fastq files.

# Log
An overview of what I have done so far. Most of this should probably be moved to the wiki, but for now I just want it in one place.

## Retrieving metadata and raw data
I searched [NCBI](https://www.ncbi.nlm.nih.gov/sra) for the accession number of the samples (PRJNA400310) from the paper. Then I downloaded the runinfo table file ‘SraRunTable.txt’ containing info on all runs and moved it to my metadata directory.
Since the files with the raw data reads are big and I'm limited on storage space, I created soft links to the PacBio and Illumina reads for scaffold 10 in the project directory from my raw_data directory instead of copying the files. The commands used for creating the soft links can be found in [retrieve_data](code/retrieve_data). 

## Canu genome assembly
I ran a Canu genome assembly on the PacBio reads in 'SRR6037732_scaffold_10.fq.gz' by running the script [genome_assembly.sh](code/genome_assembly.sh). In the script the parameter 'genomeSize' is set to 24.2m because of data from [NCBI](https://www.ncbi.nlm.nih.gov/Traces/wgs/NSDW01?display=contigs), where scaffold 10's length is stated to be 24 162 007.

The first Canu assembly job was started 03.28 on April 13th and immediately returned an error. The [output file](analysis/genome_assembly/canu/01_canu_assembly.out) states that the '-num_threads' parameter I had used was invalid, so I will use the 'maxThreads' parameter instead when re-submitting the job. A parameter that specifies the number of threads the job can use is needed in order for the job to actually use all the requested cores, since without it the job will run on one core, which is a waste of data as well as the job probably not being able to finish in time.

The second Canu job was started 06.27 on April 13th. There were only two changes made to [genome_assembly.sh](code/genome_assembly.sh); the parameter deciding how many cores the job will activley use was changed from '-num_threads' to 'maxThreads' (however still set to 4) and the prefix as well as the output directory changed name from '01_canu_assembly" to '02_canu_assembly' in order to be able to distiguish between the two jobs.
I was however nervous about accidentaly either requesting too many cores or not using all the requested cores, so I generated a graph over CPU usage for the job using the jobinfo command (full command with parameters can be found in [misc.txt](code/misc.txt)). The first graph of the job was made after 1h, 43min and 23s and can be found in [job_info_canu_01.png](job_info/job_info_canu_01.png) together with a short explanatory [text file](job_info/job_info_canu.txt). The blue graph (core usage) is consistently close to the max capacity of 400% core usage, which means that the job is almost fully using all 4 requested cores and that I am (probably and hopefully) not wasting data on this job.

The second Canu job took about 14h and 30min to run. I realised that the output files of course should be stored in the data directory and not in [analysis](analysis/), so I moved all output files there besides the .out files. The .out files contain everything that would have been printed if the script had been run in a terminal. I also removed the sub directories in [canu](analysis/genome_assembly/canu/) since each job now only will store one file here.

## Quality control and pre-processing of Illumina reads
I ran a quality control on the Illumina reads with FastQC from [fastqc_illumina.sh](code/fastqc_illumina.sh). I wasn't sure how long time it would take or how many cores I need since there were no estimations for paper V in the student manual. I assumed that my data would take as much time and use as much CPU as data from the other papers since it's the same analysis, so judging from the other papers' expected runtime for FastQC (3 min for paper II, 15 min for paper IV) I assumed it would take less than one hour for my data as well, so I submitted the job for max one hour on one core.
Running FastQC on all Illumina reads took about two minutes using one core and each Illumina read resulted in two files, one html file containing the report generated by the analysis and one zip-compressed directory that I unzipped using the unzip command, see [misc.txt](code/misc.txt).

Each unzipped directory contains two .txt files, one .fo file, one html file containing a report of the read as well as two directories with icons and images that I assume are unnecessary unless you use the graphical version of the module.
After I unzipped the directories I moved the original .zip files to a separate directory [zip_files](analysis/pre_processing/fastqc/zip_files/), although I don't think I need them any more so I could as easily have deleted them. Currently I think I'm fine on storage so I'm gonna let them be.

To summarize the output from the FastaQC analysis I used the MultiQC module to create a folder as well as an html file with data summary from the output files in [fastqc](analysis/pre_processing/fastqc/). The commands for running MultiQC can be found in [misc.txt](code/misc.txt). Illumina reads with 1P in the name are paired forward read, 1U are unpaired forward read (only when the reverse reads fails the quality control), 2P are paired reverse reads and 2U are unpaired reverse reads (only when the forward reads fails the quality control). As can be seen in [multiqc_report.html](analysis/pre_processing/fastqc/multiqc_report.html), sequences were only found in 1P and 2P, ie only in paired reads.
I don't really know what conclusion to draw from that. Does that mean that no read failed their quality control since thats the only way, at least how I understood it, that reads end up as unpaired?

If you look at the first table in the [fastqc_report](analysis/pre_processing/fastqc/SRR6058604_scaffold_10.1U_fastqc/fastqc_report.html) for 1U reads it says that the file typ analyzed was 'null' which is weird. Almost all other values are set to 0, meaning that the file contains no sequences? The only information to be found in the report is the original file name (SRR6058604_scaffold_10.1U.fastq.gz) and that it's an Illumina read.
Compared to the [fastqc_report](analysis/pre_processing/fastqc/SRR6058604_scaffold_10.1P_fastqc/fastqc_report.html) for 1P reads where all information is filled in, as I imagined all four files would have looked like.

Another thing I did to get a better overview of the results was to save all summary.txt files from the FastQC analysis in one file called [summaries_fastqc.txt](analysis/pre_processing/fastqc/summaries_fastqc.txt), see [misc.txt](code/misc.txt) for the command used.


# To-do list
* Update project plan according to feedback
* Write the wiki - or is it ok to have the same info in this file but maybe thats more difficult to organize
* Draw conclusion from the results of FastQC with the help of the MultiQC summary
* Pre-processing (trimming) of Illumina reads using Trimmomatic
* Make data_organisation.png not look like shit
