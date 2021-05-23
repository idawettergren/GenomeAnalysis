#!/bin/bash -l
#SBATCH -A g2021012
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 06:00:00
#SBATCH -J braker
#SBATCH --mail-type=ALL
#SBATCH --mail-user ida.wettergren.8542@student.uu.se 

# Load modules
module load bioinfo-tools
module load perl/5.24.1
module load GeneMark/4.33-es_Perl5.24.1
module load augustus/3.2.3_Perl5.24.1
module load bamtools/2.5.1
module load blast/2.9.0+
module load samtools/1.8
module load GenomeThreader/1.7.0
module load braker/2.1.1_Perl5.24.1

# Flags
export AUGUSTUS_CONFIG_PATH=/home/ida/genome_analysis/project/augustus/config/
export AUGUSTUS_BIN_PATH=/sw/bioinfo/augustus/3.4.0/snowy/bin/
export AUGUSTUS_SCRIPTS_PATH=/sw/bioinfo/augustus/3.4.0/snowy/scripts/
export GENEMARK_PATH=/sw/bioinfo/GeneMark/4.33-es/snowy/

# Paths
genome="/home/ida/genome_analysis/project/data/assembly/repeatmasker/improved_pacbio_assembly.fasta.masked"
aligned="/proj/g2021012/nobackup/work/ida/Aligned.sortedByCoord.out.bam"

# Move to output directory
cd /proj/g2021012/nobackup/work/ida/braker/

# Run BRAKER
braker.pl --cores=2 --softmasking --species=Durian --useexisting --genome=$genome --bam=$aligned \
--AUGUSTUS_CONFIG_PATH=/home/ida/genome_analysis/project/augustus/config/ \
--AUGUSTUS_BIN_PATH=/sw/bioinfo/augustus/3.4.0/snowy/bin/ \
--AUGUSTUS_SCRIPTS_PATH=/sw/bioinfo/augustus/3.4.0/snowy/scripts/ \
--GENEMARK_PATH=/sw/bioinfo/GeneMark/4.33-es/snowy/
