#!/bin/bash
#SBATCH --job-name=cluster_checkv

source deactivate
conda activate tools

# 1. All-vs-All BLASTN and Clustering (MIUViG parameters)
makeblastdb -in DWTP_virus_nr.fasta -dbtype nucl -out 01Clustergenome/DWTP_virus_nr
blastn -query DWTP_virus_nr.fasta -db 01Clustergenome/DWTP_virus_nr -outfmt '6 std qlen slen' -max_target_seqs 10000 -out 01Clustergenome/DWTP_virus_nr.tsv -num_threads 32 

python 01Clustergenome/anicalc.py -i 01Clustergenome/DWTP_virus_nr.tsv -o 01Clustergenome/DWTP_virus_ani.tsv
python 01Clustergenome/aniclust.py --fna DWTP_virus_nr.fasta --ani 01Clustergenome/DWTP_virus_ani.tsv --out 01Clustergenome/DWTP_virus_clusters.tsv --min_ani 95 --min_tcov 85 --min_qcov 0

cut -f 1 01Clustergenome/DWTP_virus_clusters.tsv > 01Clustergenome/DWTP_virus_clusters.txt
seqtk subseq DWTP_virus_nr.fasta 01Clustergenome/DWTP_virus_clusters.txt > 01Clustergenome/DWTP_VCs.fasta

# 2. CheckV Quality Control
conda activate vs2
checkv end_to_end 01Clustergenome/DWTP_VCs.fasta 03checkv -t 50 -d /path/to/checkv-db
# Note: Further filtering strictly retains contigs with viral genes > 0, or length > 10kb if host genes <= 1.