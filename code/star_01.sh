#!/bin/bash -l
#SBATCH -A g2021012
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 00:20:00
#SBATCH -J star
#SBATCH --mail-type=ALL
#SBATCH --mail-user ida.wettergren.8542@student.uu.se

# Load modules
module load bioinfo-tools
module load star

# Paths
genomeDir="/home/ida/genome_analysis/project/data/star/genome_index/"
masked_genome="/home/ida/genome_analysis/project/data/assembly/repeatmasker/improved_pacbio_assembly.fasta.masked"
reads_92_1="/home/ida/genome_analysis/project/data/raw_data/transcriptome/trimmed/SRR6040092_scaffold_10.1.fastq.gz"
reads_92_2="/home/ida/genome_analysis/project/data/raw_data/transcriptome/trimmed/SRR6040092_scaffold_10.2.fastq.gz"
reads_93_1="/home/ida/genome_analysis/project/data/raw_data/transcriptome/trimmed/SRR6040093_scaffold_10.1.fastq.gz"
reads_93_2="/home/ida/genome_analysis/project/data/raw_data/transcriptome/trimmed/SRR6040093_scaffold_10.2.fastq.gz"
reads_94_1="/home/ida/genome_analysis/project/data/raw_data/transcriptome/trimmed/SRR6040094_scaffold_10.1.fastq.gz"
reads_94_2="/home/ida/genome_analysis/project/data/raw_data/transcriptome/trimmed/SRR6040094_scaffold_10.2.fastq.gz"
reads_95_1="/home/ida/genome_analysis/project/data/trimmomatic/SRR6040095_scaffold_10.1P.fastq.gz"
reads_95_2="/home/ida/genome_analysis/project/data/trimmomatic/SRR6040095_scaffold_10.2P.fastq.gz"
reads_96_1="/home/ida/genome_analysis/project/data/raw_data/transcriptome/trimmed/SRR6040096_scaffold_10.1.fastq.gz"
reads_96_2="/home/ida/genome_analysis/project/data/raw_data/transcriptome/trimmed/SRR6040096_scaffold_10.2.fastq.gz"
reads_97_1="/home/ida/genome_analysis/project/data/raw_data/transcriptome/trimmed/SRR6040097_scaffold_10.1.fastq.gz"
reads_97_2="/home/ida/genome_analysis/project/data/raw_data/transcriptome/trimmed/SRR6040097_scaffold_10.2.fastq.gz"
reads_66_1="/home/ida/genome_analysis/project/data/raw_data/transcriptome/trimmed/SRR6156066_scaffold_10.1.fastq.gz"
reads_66_2="/home/ida/genome_analysis/project/data/raw_data/transcriptome/trimmed/SRR6156066_scaffold_10.2.fastq.gz"
reads_67_1="/home/ida/genome_analysis/project/data/raw_data/transcriptome/trimmed/SRR6156067_scaffold_10.1.fastq.gz"
reads_67_2="/home/ida/genome_analysis/project/data/raw_data/transcriptome/trimmed/SRR6156067_scaffold_10.2.fastq.gz"
reads_69_1="/home/ida/genome_analysis/project/data/raw_data/transcriptome/trimmed/SRR6156069_scaffold_10.1.fastq.gz"
reads_69_2="/home/ida/genome_analysis/project/data/raw_data/transcriptome/trimmed/SRR6156069_scaffold_10.2.fastq.gz"

# Move to output directory
cd /proj/g2021012/nobackup/work/ida/

# Create index with STAR
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $masked_genome --genomeSAindexNbases 9

# Map RNA reads to masked assembly with STAR
STAR --runThreadN 8 --runMode alignReads --genomeDir $genomeDir --readFilesCommand gunzip -c --readFilesIn $reads_92_1,$reads_92_2,\
$reads_93_1,$reads_93_2,$reads_94_1,$reads_94_2,$reads_95_1,$reads_95_2,$reads_96_1,$reads_96_2,$reads_97_1,$reads_97_2,$reads_66_1,\
$reads_66_2,$reads_67_1,$reads_67_2,$reads_69_1,$reads_69_2 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 4010136799
