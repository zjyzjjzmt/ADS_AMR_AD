#!/bin/bash
#SBATCH --job-name=functional_profiling
#SBATCH --cpus-per-task=50

source deactivate

# 1. Lifestyle Prediction: VIBRANT
echo "[INFO] Running VIBRANT for lifestyle prediction..."
conda activate tools
VIBRANT_run.py -i 01Clustergenome/DWTP_VCs.fasta \
               -f nucl -t 50 -virome \
               -folder vibrant_out

# 2. Lifestyle Prediction: PhaTYP (PhaBOX2)
echo "[INFO] Running PhaTYP for lifestyle prediction..."
conda deactivate
conda activate phabox2
phabox2 --task phatyp \
        --dbdir /path/to/phabox_db \
        --outpth phatyp_out \
        --contigs 01Clustergenome/DWTP_VCs.fasta \
        --threads 50

# 3. ORF Prediction and ARG/VF Mapping 
echo "[INFO] Predicting ORFs and mapping to functional databases..."
conda deactivate
conda activate tools

# Predict ORFs
prodigal -i 01Clustergenome/DWTP_VCs.fasta \
         -o DWTP_VCs.genes \
         -a DWTP_VCs.faa \
         -d DWTP_VCs.fna -p meta

mkdir -p 06ARGs 07VFs

# ARG Mapping (Strict parameters: evalue 1e-10, cov 70%, id 80%)
diamond blastx -d path/to/SARG.2.2_nr -q DWTP_VCs.fna -o 06ARGs/DWTP_VCs-SARG.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
diamond blastx -d path/to/card3.1.4_nr -q DWTP_VCs.fna -o 06ARGs/DWTP_VCs-card.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1

# Virulence Factor (VF) Mapping
diamond blastx -d path/to/VFDB_nr -q DWTP_VCs.fna -o 07VFs/DWTP_VCs-VFDB.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1

echo "[INFO] Functional profiling completed."