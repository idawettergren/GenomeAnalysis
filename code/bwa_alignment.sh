#!/bin/bash -l
#SBATCH -A g2021012
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J bwa_alignment
#SBATCH --mail-type=ALL
#SBATCH --mail-user ida.wettergren.8542@student.uu.se

# Load modules
module load bioinfo-tools
module load bwa/0.17.7

# Paths to the data
illumina_1="/home/ida/genome_analysis/project/data/raw_data/illumina/SRR6058604_scaffold_10.1P.fastq.gz"
illumina_2="/home/ida/genome_analysis/project/data/raw_data/illumina/SRR6058604_scaffold_10.2P.fastq.gz"
assembly="/home/ida/genome_analysis/project/data/assemblies/pacbio_canu_assembly/02_canu_assembly.contigs.fasta"

# Move to output directory
cd /home/ida/genome_analysis/project/GenomeAnalysis/analysis/alignment/

# Run BWA
