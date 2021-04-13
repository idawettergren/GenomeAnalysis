# I think this is fine to run on the login node ??

# load modules
module load bioinfo-tools
module load MultiQC/1.9

# summarize FastQC results
multiqc /home/ida/genome_analysis/project/GenomeAnalysis/analysis/pre_processing/fastqc/

