#!/bin/bash
#SBATCH --job-name=ont_viral_discovery

source deactivate
conda activate flye

# 1. Flye Assembly (Meta mode)
for sample in 01RAW_READS/*.fastq; do
    sample_name=$(basename $sample .fastq)
    flye --nano-raw $sample --out-dir Flye_Assemblies/$sample_name --genome-size 5g --meta --threads 40
done

# 2. Viral Prediction via geNomad
conda activate genomad
for assembly in Flye_Assemblies/*/assembly.fasta; do
    sample_name=$(basename $(dirname $assembly))
    genomad end-to-end --cleanup $assembly 03Virus_Prediction/$sample_name /path/to/genomad_db
done

cat 03Virus_Prediction/*/assembly_summary/assembly_virus.fna > 03Virus_Prediction/all_predicted_viruses.fasta

# 3. Merge Illumina and ONT Viral Candidates
conda activate quickmerge
merge_wrapper.py 03Virus_Prediction/all_predicted_viruses.fasta 03_Viral_Discovery/illumina_viral_nr.fasta