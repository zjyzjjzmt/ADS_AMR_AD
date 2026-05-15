#!/bin/bash

source deactivate
conda activate tools

# Make directories
mkdir -p 04reads_based_analysis/02card
cd 04reads_based_analysis/02card

# 1. Build DIAMOND database for CARD
echo "[INFO] Building DIAMOND database..."
diamond makedb --in CARD3.1.4.fasta -d card3.1.4_nr

cd ../../
# 2. Alignment using DIAMOND
echo "[INFO] Starting DIAMOND blastx alignment..."
for i in 01Reads_fq/*.fastq.gz; do
    diamond blastx -d 04reads_based_analysis/02card/card3.1.4_nr \
                   -q "$i" -o "$i-CARD.txt" \
                   --evalue 1e-5 --query-cover 75 --id 90 -k 1 --threads 60
    mv "$i-CARD.txt" 04reads_based_analysis/02card/
done

# 3. Generating raw count table using Awk
cd 04reads_based_analysis/02card
echo "[INFO] Compiling hit counts..."
awk -F'\t' '{print $1"|"$2"|"$3"|"$4}' CARD3.1.4.txt > id_list.txt

# Create headers and initialize output
FILES=(*.fastq.gz-CARD.txt)
{
    echo -ne "ID"
    for f in "${FILES[@]}"; do echo -ne "\t$(basename "$f")"; done
    echo
} > result.tsv

# Index and sum counts in memory
ls *.fastq.gz-CARD.txt > file_list.tmp
awk -v idfile="id_list.txt" -v filelist="file_list.tmp" '
BEGIN {
    while ((getline < idfile) > 0) target[$1] = 1
    file_count = 0
    while ((getline f < filelist) > 0) files[++file_count] = f
    
    for (i = 1; i <= file_count; i++) {
        file = files[i]
        while ((getline < file) > 0) {
            id = $2
            if (id in target) count[file SUBSEP id]++
        }
        close(file)
    }
    
    for (id in target) {
        line = id
        for (i = 1; i <= file_count; i++) {
            c = count[files[i] SUBSEP id] + 0
            line = line "\t" c
        }
        print line
    }
}' >> result.tsv
rm -f file_list.tmp

# Remove rows with only zeros
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {print} NR>1 {sum=0; for (i=2; i<=NF; i++) sum+=$i; if (sum!=0) print}' result.tsv > count_table_final.tsv

# Merge with CARD metadata
csvtk join -t --left-join -f 1 count_table_final.tsv CARD_mapping-20220905.txt > Study_CARD_mapping.tsv
echo "[INFO] ARG mapping completely finished."