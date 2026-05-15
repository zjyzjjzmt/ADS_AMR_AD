#!/bin/bash
#SBATCH --job-name=rna_clustering
#SBATCH --cpus-per-task=50

source deactivate
conda activate tools
cd 08RNA_VCs
mkdir -p 01Clustergenome 03checkv

# 1. All-vs-all BLAST and Clustering (90% ANI, 80% Coverage)
echo "[INFO] Clustering RNA vOTUs..."
makeblastdb -in AD_RNA_virus_nr.fasta -dbtype nucl -out 01Clustergenome/AD_RNA_virus_nr
blastn -query AD_RNA_virus_nr.fasta -db 01Clustergenome/AD_RNA_virus_nr -outfmt '6 std qlen slen' -max_target_seqs 10000 -out 01Clustergenome/AD_RNA_virus_nr.tsv -num_threads 32 

python 01Clustergenome/anicalc.py -i 01Clustergenome/AD_RNA_virus_nr.tsv -o 01Clustergenome/AD_RNA_virus_ani.tsv
python 01Clustergenome/aniclust.py --fna AD_RNA_virus_nr.fasta --ani 01Clustergenome/AD_RNA_virus_ani.tsv --out 01Clustergenome/AD_RNA_virus_clusters.tsv --min_ani 90 --min_tcov 80 --min_qcov 0

cut -f 1 01Clustergenome/AD_RNA_virus_clusters.tsv > 01Clustergenome/AD_RNA_virus_clusters.txt
seqtk subseq AD_RNA_virus_nr.fasta 01Clustergenome/AD_RNA_virus_clusters.txt > 01Clustergenome/AD_RNA_VCs.fasta

# 2. CheckV Quality Control
echo "[INFO] Running CheckV QC..."
conda deactivate
conda activate vs2
checkv end_to_end 01Clustergenome/AD_RNA_VCs.fasta 03checkv -t 50 -d /path/to/checkv-db