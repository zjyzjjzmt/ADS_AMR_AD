#!/bin/bash
#SBATCH --job-name=rna_annotation

# 1. ARG Mapping
conda activate tools
prodigal -i 01Clustergenome/AD_RNA_VCs.fasta -o 01Clustergenome/AD_RNA_VCs.genes -a 01Clustergenome/AD_RNA_VCs.faa -d 01Clustergenome/AD_RNA_VCs.fna -p meta

mkdir -p 06ARGs
diamond blastx -d /path/to/SARG.2.2_nr -q 01Clustergenome/AD_RNA_VCs.fna -o 06ARGs/AD_RNA_VCs-SARG.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
diamond blastx -d /path/to/card3.1.4_nr -q 01Clustergenome/AD_RNA_VCs.fna -o 06ARGs/AD_RNA_VCs-card.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1

# 2. Taxonomy (geNomad, CAT, vConTACT2)
# Executed standard classification pipelines identically to the DNA virus module.
# (Code abridged for brevity; utilizes geNomad end-to-end, CAT contigs, and vcontact2_gene2genome)