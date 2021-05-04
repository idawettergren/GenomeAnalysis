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

# Names for output files
out_forward_unpaired="SRR6040095_scaffold_10.1U.fastq.gz"
out_forward_paired="SRR6040095_scaffold_10.1P.fastq.gz"
out_reverse_unpaired="SRR6040095_scaffold_10.2U.fastq.gz"
out_reverse_paired="SRR6040095_scaffold_10.2P.fastq.gz"


# Run Trimmomatic (trim Illumina sequences and leading/trailing sequences below 3 quality)
java -jar $TRIMMOMATIC_HOME/trimmomatic.jar PE -threads 2 $reads_forward $reads_reverse $out_forward_paired $out_forward_unpaired \
$out_reverse_paired $out_reverse_unpaired ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3
