#!/bin/bash
#SBATCH --job-name=rna_host_prediction
#SBATCH --cpus-per-task=50

source deactivate
conda activate tools
mkdir -p 10host-link/01crispr 10host-link/02EVEs 10host-link/03RNAVirHost

# 1. EVEs (Endogenous Viral Elements) Linkage via tBLASTn
echo "[INFO] Running EVEs mapping..."
makeblastdb -dbtype nucl -in /path/to/MAGs_database.fa -input_type fasta -title MAGs_DB -out 10host-link/02EVEs/MAGs_DB -blastdb_version 5
tblastn -query 10host-link/02EVEs/RDRP.fa -db 10host-link/02EVEs/MAGs_DB -out 10host-link/02EVEs/MAGs_hosts_EVEs.txt \
        -evalue 1e-20 -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident bitscore mismatch qcovs qcovhsp" -num_threads 50

# 2. CRISPR-Cas Spacer Matches (blastn-short)
echo "[INFO] Running CRISPR spacer matching..."
makeblastdb -dbtype nucl -in 01Clustergenome/AD_RNA_VCs.fasta -input_type fasta -title AD_RNA_VCs -out 10host-link/01crispr/AD_RNA_VCs -blastdb_version 5
# Stringent parameters for short spacer alignments: 90% ID, word_size 7, dust off
blastn -task blastn-short -query /path/to/MAGs_spacers_nr.fasta -db 10host-link/01crispr/AD_RNA_VCs \
       -perc_identity 90 -out 10host-link/01crispr/RNA_VCs_MAGs_criprcas_host.txt \
       -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident bitscore mismatch qcovs qcovhsp" \
       -max_target_seqs 1 -dust no -word_size 7 -num_threads 50

# 3. ML-based Taxonomy Prediction using RNAVirHost
echo "[INFO] Running RNAVirHost prediction..."
conda deactivate
conda activate RNAVirHost

rnavirhost classify_order -i 01Clustergenome/AD_RNA_VCs.fasta -o 10host-link/03RNAVirHost/RNAVirHost_TAXONOMIC
rnavirhost predict -i 01Clustergenome/AD_RNA_VCs.fasta --taxa 10host-link/03RNAVirHost/RNAVirHost_TAXONOMIC -o 10host-link/03RNAVirHost/RNAVirHost