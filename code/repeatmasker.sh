#!/bin/bash -l
#SBATCH -A g2021012
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J repeatmasker
#SBATCH --mail-type=ALL
#SBATCH --mail-user ida.wettergren.8542@student.uu.se

# Load modules
module load bioinfo-tools
module load RepeatMasker

# Paths
outdir="/home/ida/genome_analysis/project/GenomeAnalysis/analysis/repeatmasker/"
assembly="/home/ida/genome_analysis/project/data/assembly/pilon_improvement/improved_pacbio_assembly.fasta"

# Run RepeatMasker
RepeatMasker -dir $outdir -xsmall $assembly
