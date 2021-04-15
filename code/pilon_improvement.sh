#!/bin/bash -l
#SBATCH -A g2021012
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J pilon_improvement
#SBATCH --mail-type=ALL
#SBATCH --mail-user ida.wettergren.8542@student.uu.se

# Load modules
module load bioinfo-tools
module load samtools
module load Pilon/1.24

# Move to directory with .bam file
cd /home/ida/genome_analysis/project/data/alignments/

# Sort .bam file and create index
samtools sort -@ 2 illumina_pacbio_alignment.bam -o illumina_pacbio_alignment_sorted.bam
samtools index illumina_pacbio_alignment_sorted.bam

# Paths to files
assembly="/home/ida/genome_analysis/project/data/assemblies/pacbio_canu_assembly/02_canu_assembly.contigs.fasta"
sorted_bam="/home/ida/genome_analysis/project/data/alignments/illumina_pacbio_alignment_sorted.bam"

# Move to output directory
cd /home/ida/genome_analysis/project/data/assemblies/pilon_improvement/

# Run Pilon on the sorted assembly to improve it


