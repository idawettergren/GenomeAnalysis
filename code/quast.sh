#!/bin/bash -l
#SBATCH -A g2021012
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J quast_quality_control
#SBATCH --mail-type=ALL
#SBATCH --mail-user ida.wettergren.8542@student.uu.se

# Load modules
module load bioinfo-tools
module load quast

# Paths to files (assembly after pilon improvements, reference genome)
assembly="/home/ida/genome_analysis/project/data/assembly/pilon_improvement/improved_pacbio_assembly.fasta"
out_dir="/home/ida/genome_analysis/project/GenomeAnalysis/analysis/quast/"

quast.py -t 2 -o $out_dir $assembly
