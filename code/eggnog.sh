#!/bin/bash -l
#SBATCH -A uppmaxg2022-2-5
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

emapper.py -i ?? --itype genome -o func_annotation --output_dir ?? --cpu 2
