#!/bin/bash -l
#SBATCH -A uppmax2022-2-5
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 08:00:00
#SBATCH -J eggnog
#SBATCH --mail-type=ALL
#SBATCH --mail-user ida.wettergren.8542@student.uu.se 

# Load modules
module load bioinfo-tools
module load eggNOG-mapper

emapper.py -i /home/ida/genome_analysis/project/data/assembly/repeatmasker/improved_pacbio_assembly.fasta.masked --itype genome -o func_annotation --output_dir /home/ida/genome_analysis/project/data/eggnog --cpu 2
