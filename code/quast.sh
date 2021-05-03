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
module bioinfo-tools
module load quast/5.0.2

# Paths to files (assembly after pilon improvements, reference genome)
assembly="/home/ida/genome_analysis/project/data/assembly/pilon_improvement/improved_pacbio_assembly.fasta"
ref=""
out_dir="/home/ida/genome_analysis/project/GenomeAnalysis/analysis/quast/"

quast.py -o $out_dir -r $ref $assembly
