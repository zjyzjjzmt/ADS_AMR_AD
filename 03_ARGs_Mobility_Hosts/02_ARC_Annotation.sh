#!/bin/bash
#SBATCH --job-name=ARC_annotation

source deactivate
conda activate tools

cd 07ARCs_analysis

# 1. Remove redundancy using CD-HIT-EST (100% identity, alignment coverage control)
echo "[INFO] Removing redundancy of ARCs..."
cd-hit-est -i DWTP_ARCs.fasta -o DWTP_ARCs_nr.fasta -c 1.0 -n 10 -M 0 -T 8

# Repredict ORFs on the non-redundant catalog
prodigal -i DWTP_ARCs_nr.fasta -o DWTP_ARCs_nr.genes -a DWTP_ARCs_nr.pro.fa -d DWTP_ARCs_nr.nul.fa -p meta

# 2. Comprehensive Annotation (CARD, SARG, ICEs, MGEs)
echo "[INFO] Annotating non-redundant ARCs..."
mkdir -p 01SARG 02card 04ICEs 05MGEs

# ARG (SARG & CARD)
diamond blastx -d path/to/SARG.2.2_nr -q DWTP_ARCs_nr.nul.fa -o 01SARG/DWTP_ARCs-SARG.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
cut -f 1,2 01SARG/DWTP_ARCs-SARG.txt | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' | cut -f 1,2,3,4,5,6,8 | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' | sort -k 1n | datamash -sW -g1 collapse 2 > 01SARG/DWTP_ARCs-SARG_list2.txt

diamond blastx -d path/to/card3.1.4_nr -q DWTP_ARCs_nr.nul.fa -o 02card/DWTP_ARCs-card.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
cut -f 1,2 02card/DWTP_ARCs-card.txt | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' | cut -f 1,2,3,4,5,6,8 | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' | sort -k 1n | awk '{split ($2, T, "|"); $2 = T[4]}1' OFS="\t" | datamash -sW -g1 collapse 2 > 02card/DWTP_ARCs-card_list2.txt

# ICEs (Integrative and Conjugative Elements)
diamond blastx -d path/to/ICE3.0_nr -q DWTP_ARCs_nr.nul.fa -o 04ICEs/DWTP_ARCs-ICE.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
cut -f 1,2 04ICEs/DWTP_ARCs-ICE.txt | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' | cut -f 1,2,3,4,5,6,8 | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' | sort -k 1n | datamash -sW -g1 collapse 2 > 04ICEs/DWTP_ARCs-ICE_list2.txt

# MGEs (Mobile Genetic Elements via BLASTN)
blastn -query DWTP_ARCs_nr.nul.fa -db path/to/MGEs -out 05MGEs/DWTP_ARCs-MGEs.txt -evalue 1e-10 -perc_identity 80 -num_threads 10 -outfmt 6 -max_target_seqs 1
cut -f 1,2 05MGEs/DWTP_ARCs-MGEs.txt | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' | cut -f 1,2,3,4,5,6,8 | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' | sort -k 1n | datamash -sW -g1 collapse 2 > 05MGEs/DWTP_ARCs-MGEs_list2.txt

echo "[INFO] Annotation finished."