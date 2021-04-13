#!/bin/bash -l

#SBATCH -A g2021012
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J fastqc
#SBATCH --mail-type=ALL
#SBATCH --mail-user ida.wettergren.8542@student.uu.se

# Load FastQC module
module load bioinfo-tools
module load FastQC/0.11.9

# Run FastQC on all Illumina reads in the directory
fastqc -o /home/ida/genome_analysis/project/GenomeAnalysis/analysis/pre_processing/01_fastqc/ \
/home/ida/genome_analysis/project/data/raw_data/illumina/*gz



