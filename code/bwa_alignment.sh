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
module load bwa/0.7.17
module load samtools

# Paths to the data
illumina_1="/home/ida/genome_analysis/project/data/raw_data/illumina/SRR6058604_scaffold_10.1P.fastq.gz"
illumina_2="/home/ida/genome_analysis/project/data/raw_data/illumina/SRR6058604_scaffold_10.2P.fastq.gz"
assembly="/home/ida/genome_analysis/project/data/assemblies/pacbio_canu/02_canu_assembly.contigs.fasta"

# Move to output directory
cd /home/ida/genome_analysis/project/data/alignments/illumina_pacbio/

# Align BWA with Illumina short reads
# Create index
bwa index $assembly

# Align each Illumina read
bwa aln -t 2 $assembly $illumina_1 > aln_illumina_1.sai
bwa aln -t 2 $assembly $illumina_2 > aln_illumina_2.sai

# Generate alignment from the (paired-end) reads
bwa sampe $assembly aln_illumina_1.sai aln_illumina_2.sai $illumina_1 $illumina_2 > illumina_pacbio_alignment.sam

# Convert .sam file to .bam
samtools view -b illumina_pacbio_alignment.sam > illumina_pacbio_alignment.bam
