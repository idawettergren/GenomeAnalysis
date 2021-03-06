#!/bin/bash -l

#SBATCH -A g2021012
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 20:00:00
#SBATCH -J canu_assembly
#SBATCH --mail-type=ALL
#SBATCH --mail-user ida.wettergren.8542@student.uu.se

# Load Canu module
module load bioinfo-tools
module load canu/2.0

# Run Canu
canu -maxThreads=4 -p 02_canu_assembly -d \
/home/ida/genome_analysis/project/GenomeAnalysis/analysis/assembly/pacbio_canu/02_canu_assembly/ \
genomeSize=24.2m -pacbio-raw \
/home/ida/genome_analysis/project/data/raw_data/pacbio/SRR6037732_scaffold_10.fq.gz

