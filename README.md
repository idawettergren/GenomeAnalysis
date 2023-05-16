# Genome Analysis

## Aim of project
The aim of this project is to reproduce some of the results from ["The draft genome of tropical fruit durian (Durio zibethinus)" by Bin Tean Teh et al](https://www.nature.com/articles/ng.3972/). Scaffold 10 will be analysed instead of the entire genome.

## Methods
To assemble the provided PacBio reads Canu will be used. For quality control of provided Illumina reads FastQC will be used, after which Trimmomatic will be used to trim the reads if necessary. The Illumina reads will then be aligned with the PacBio assembly using BWA and can then be used to improve the assembly using Pilon.

To find where in the assembly potential genes are BRAKER will be used to perform a structural annotation, and to find potential functions of the genes eggNOGmapper will be used to perform a functional annotation.

Finally Htseq and Deseq2 will be used to determine up- and down regulation of the genes.

## Timeline
In order to be on track the genome assembly should be done on April 18th, the annotation on May 2nd and the RNA mapping on May 9th. The timeline from the student manual have been moved up slightly here so that each deadline is on a Sunday, which allows for work on the weekends, should it be needed.

## Data organisation
It's unclear how much space will be needed so the goal will be to keep files compressed where possible and to use soft links for raw data to conserve space. All files used in the project will be stored in my UPPMAX directory. This git repository will contain all files except for those that contain data. These files will be many different formats, such as .fasta, .fastq or .bam.

The files [data_organisation.txt](data_organisation.txt) and [data_organisation.png](data_organisation.png) give an overview of the UPPMAX directory used to store the files.

# Pipeline
1. Run pacbio_canu_assembly.sh to assemble the PacBio reads with Canu. Output: /analysis/assembly/pacbio_canu/02_canu_assembly/
2. Run fastqc_illumina.sh to perform a quality control with FastQC on the Illumina reads. Output: /analysis/pre_processing/fastqc/illumina/
3. Run bwa_alignment.sh to align the PacBio assembly with the Illumina reads using BWA. Output: /data/alignment/illumina_pacbio/
4. Run pilon_improvement.sh to improve the assembly using the aligned Illumina reads with Pilon. Output: /data/assembly/pilon_improvement/
5. Run quast.sh to evalute the Pilon improvement. Output: /analysis/quast/
6. Run trimmomatic.sh to trim the RNA reads. Output: /data/raw_data/transcriptome/trimmed/SRR6040095_scaffold_10.1U.fastq.gz (.1P/2U/2P)
7. Run fastqc_rna.sh to perfrom a quality control with FastQC on the RNA reads. Output: /analysis/pre_processing/fastqc/rna/
8. Run star_01.sh to map the RNA reads to the Pilon improved assembly. Output: /proj/g2021012/nobackup/work/ida/
9. Run braker.sh.
10. Run EggNogMapper. (Currently here)
11. Run HTSeq.
12. Run Deseq2.

# Log
An overview of what I have done and my thoughts surrounding the project.

## Retrieving data
I searched [NCBI](https://www.ncbi.nlm.nih.gov/sra) for the accession number of the samples (PRJNA400310) from the paper and downloaded the runinfo table file ‘SraRunTable.txt’ containing info on all runs and moved it to my metadata directory.

Since the provided files containing the raw data are large and I'm limited on storage space, I created soft links to the PacBio, Illumina and RNA reads (for scaffold 10) in the project directory from my raw_data directory instead of copying the files there. The commands used for creating the soft links can be found in [misc.txt](code/misc.txt). 

## Quality control and pre-processing of Illumina reads
I ran a quality control on the provided Illumina reads using FastQC in the script [fastqc_illumina.sh](code/fastqc_illumina.sh). I wasn't sure how long time it would take or how many cores I need since there were no estimations for paper V in the student manual, so I assumed that it would be similar to the other papers'. Paper II's expected runtime was 3 min and paper IV's was 15 min, so I submitted the job for max one hour using one core.

Running FastQC took about two minutes and each of the four reads tested returned two files; one .html file containing the report of the analysis and one zip-compressed directory that I unzipped using the unzip command, see [misc.txt](code/misc.txt). It might have been a mistake to unzip the compressed directory since I think MultiQC can read compressed directories.

Each unzipped directory contains two .txt files, one .fo file, one .html file containing a report of the tests as well as two directories with icons and images that I assume are unnecessary unless you use the graphical version of the module.
After this the .zip files were moved to a separate directory [zip_files](analysis/pre_processing/fastqc/illumina/zip_files/).

To summarize the FastaQC results for the four reads in one file and get an overview of the tests, I ran MultiQC on the output files and saved all summary.txt files in one file called [summaries_fastqc.txt](analysis/pre_processing/fastqc/illumina/summaries_fastqc.txt), see [misc.txt](code/misc.txt) for the commands used. The resulting report from MultiQC can be found in [multiqc_report.html](analysis/pre_processing/fastqc/illumina/multiqc/multiqc_report.html). Illumina reads with 1P in the name are paired forward reads, 1U are unpaired forward reads, 2P are paired reverse reads and 2U are unpaired reverse reads.

The provided paired Illumina reads appear to already be pre-processed, since no adapters were found and the sequences had high quality scores. I will therefore not trim them with Trimmomatic. I could run Trimmomatic followed by another FastQC analysis and compare the results to the FastQC analysis in order to see if there was any significant improvements gained from trimming the reads. However, since they already passed the quality control I will use the provided reads directly to correct the PacBio assembly.

## Assembly of PacBio reads
I ran a Canu genome assembly on the PacBio reads in 'SRR6037732_scaffold_10.fq.gz' by running the script [genome_assembly.sh](code/pacbio_canu_assembly.sh). The parameter 'genomeSize' is set to 24.2m because of data from [NCBI](https://www.ncbi.nlm.nih.gov/Traces/wgs/NSDW01?display=contigs), where scaffold 10's length is stated to be 24 162 007.

In the student manual the estimated time for this job was 17 hour when using four cores, so I submitted it for max 20 hours using four cores. 

The first Canu assembly job was started on April 13th and immediately returned an error. The [.out file](analysis/assembly/pacbio_canu/01_canu_assembly.out) states that the '-num_threads' parameter used was invalid, so when re-submitting the job the 'maxThreads' parameter will be used instead. A parameter that specifies the number of threads the job can use is needed in order for the job to actually use all the requested cores. Without it the job will run on one core, which is a waste of data and will probably result in the job not being completed in time.

The second Canu job was started later the same day. There were two changes made to [pacbio_canu_assembly.sh](code/pacbio_canu_assembly.sh): the parameter allowing the job to use multiple cores was changed from '-num_threads' to 'maxThreads' (however still set to 4) and the prefix as well as the output directory changed name from '01_canu_assembly" to '02_canu_assembly' in order to be able to distiguish between outputs from the two jobs.

I was a bit paranoid about accidentally either requesting too many cores or not using the requested cores, so I generated a graph over CPU usage 1h, 43min and 23s into the job, using the jobinfo command in [misc.txt](code/misc.txt). The graph can be found in [job_info_canu_01.png](job_info/job_info_canu_01.png) together with a short explanatory [text file](job_info/job_info_canu.txt). The blue graph (core usage) is consistently close to the max capacity of 400% core usage, which means that the job is almost fully using all 4 requested cores and that I am (probably and hopefully) not wasting data on this job. 

The second Canu job took about 14h and 30min to complete. I realised that I wanted the output files to be stored in the data directory, so I moved them all (besides the two .out files). .out files contain everything that would have been printed if the script had been run in a terminal. I also removed the sub directories in [pacbio_canu](analysis/assembly/pacbio_canu/), since each job now only store one file each here.

## Alignment of PacBio assembly and Illumina reads
To be able to correct the Canu assembly of the PacBio reads using the Illumina reads, they first have to be aligned to the assembly. This was done by running BWA on the assembly and the Illumina reads in [bwa_alignment.sh](code/bwa_alignment.sh).

First the script moves to the directory for the output, which is something I could have implemented earlier in the [FastQC](code/fastqc_illumina.sh) and [Canu](code/pacbio_canu_assembly.sh) scripts so they could be run from any directory.

The .sam file that BWA returns is converted to .bam format containing the same information, but encoded in binary so the file is smaller.

In the student manual the estimated time for this job was one hour when using two cores, so I submitted it for max two hours using two cores. It took 31min for the job to complete.

## Improve PacBio assembly with Illumina reads
To improve the PacBio assembly with the aligned Illumina reads I will use Pilon. Pilon takes a .fasta file together with a sorted .bam file of aligned reads, compares them and improve sequences where the assembly and Illumina reads differ. The output will be a .fasta file containing the improved genome sequence. The script to run this job is [pilon_improvement.sh](code/pilon_improvement.sh).

Since Pilon needs the .bam file to be sorted and the .bam file converted from the .sam file that the BWA alignment returned isn't, the script first sorts the .bam file with 'samtools sort'. This chould have been done in a pipe when running the BWA alignment, perhaps like this:

bwa sampe $assembly aln_illumina_1.sai aln_illumina_2.sai $illumina_1 $illumina_2 > illumina_pacbio_alignment.sam | samtools view -b -@ 1 illumina_pacbio_alignment.sam > illumina_pacbio_alignment.bam | samtools sort -@ 2 illumina_pacbio_alignment.bam -o illumina_pacbio_alignment_sorted.bam

This would mean that the .sam output from BWA would be piped through and converted to a .bam file which would be piped through again and sorted. This would only create the final sorted .bam file and not the intermediary ones which would have been nice.

In order to track changes made by Pilon I will use the argument '--changes' to output a .change file with all changes made in addition to the .fasta file.

## Mask repeated sequences
Ran RepeatMasker.

# Quality control and pre-processing of RNA reads
Both trimmed and untrimmed RNA reads were provided, so I decided to trim the reads myself. This was done using Trimmomatic, where the Illumina tags and reads with a quality score under 3 were removed. After the trimming of the reads, I ran a quality control using FastQC, just as I did for the DNA reads.

# Mapping of RNA reads
To map the RNA reads after the pre-processing I will use STAR.

# Annotation
To annotate the assembly, I used BRAKER. However, when running BRAKER I didn't get the output files that I wanted. The .gtf hints file was deleted after running BRAKER because it was empty and the .gff hints file was also empty. Because of this I downloaded the full .gff file from NCBI (https://www.ncbi.nlm.nih.gov/assembly/GCF_002303985.1) and will extract scaffold 10 (NW_019167827.1) using the function GFF_split_into_scaffolds from the PopGenome package in R.
