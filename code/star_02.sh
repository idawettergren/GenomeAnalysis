#!/bin/bash -l
#SBATCH -A g2021012
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 00:20:00
#SBATCH -J star
#SBATCH --mail-type=ALL
#SBATCH --mail-user ida.wettergren.8542@student.uu.se

# Load modules
module load bioinfo-tools
module load star

# Paths
genomeDir="/home/ida/genome_analysis/project/data/star/genome_index2/"
ref_genome=""
reads_95_1="/home/ida/genome_analysis/project/data/trimmomatic/SRR6040095_scaffold_10.1P.fastq.gz"
reads_95_2="/home/ida/genome_analysis/project/data/trimmomatic/SRR6040095_scaffold_10.2P.fastq.gz"

# Move to output directory
cd /home/ida/genome_analysis/project/data/star/htseq/

# Create index with STAR
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $ref_genome --genomeSAindexNbases 9

# Map RNA reads to masked assembly with STAR
STAR --runThreadN 8 --runMode alignReads --genomeDir $genomeDir --readFilesCommand gunzip -c --readFilesIn $reads_95_1 $reads_95_2 \
--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 2077397127
