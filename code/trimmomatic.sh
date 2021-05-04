#!/bin/bash -l
#SBATCH -A g2021012
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J trimmomatic_rna
#SBATCH --mail-type=ALL
#SBATCH --mail-user ida.wettergren.8542@student.uu.se

# Load modules
modules load bioinfo-tools
modules load trimmomatic

# Paths to files
reads_forward="/home/ida/genome_analysis/project/data/raw_data/transcriptome/untrimmed/SRR6040095_scaffold_10.1.fastq.gz"
reads_reverse="/home/ida/genome_analysis/project/data/raw_data/transcriptome/untrimmed/SRR6040095_scaffold_10.2.fastq.gz"
