#!/bin/bash

source deactivate
conda activate metaphlan4

MPA_DB="/path/to/db/metaphlan4"
TMP_DIR="./tmp"
mkdir -p 10metaphlan4 $TMP_DIR

echo "[INFO] Running MetaPhlAn 4..."
for F in 01Reads_fq/*_1.fastq.gz; do
    R=${F%_*}_2.fastq.gz
    BASE=${F##*/}
    SAMPLE=${BASE%_*}
    
    metaphlan $F,$R --input_type fastq \
                    --nproc 60 \
                    --offline \
                    --bowtie2out 10metaphlan4/$SAMPLE.bowtie2out \
                    --output_file 10metaphlan4/$SAMPLE.profile.txt \
                    --bowtie2db $MPA_DB \
                    --tmp_dir $TMP_DIR
done

# Merge profiles into a single abundance table
merge_metaphlan_tables.py 10metaphlan4/*.profile.txt > 10metaphlan4/merged_abundance_table.tsv

# Calculate Diversity Metrics
cd 10metaphlan4
DIV_SCRIPT="$CONDA_PREFIX/lib/python3.7/site-packages/metaphlan/utils/calculate_diversity.R"

echo "[INFO] Calculating alpha and beta diversity..."
Rscript $DIV_SCRIPT -f merged_abundance_table.tsv -d beta -m bray-curtis
Rscript $DIV_SCRIPT -f merged_abundance_table.tsv -d alpha -m richness
Rscript $DIV_SCRIPT -f merged_abundance_table.tsv -d alpha -m shannon
Rscript $DIV_SCRIPT -f merged_abundance_table.tsv -d alpha -m simpson
Rscript $DIV_SCRIPT -f merged_abundance_table.tsv -d alpha -m gini

echo "[INFO] Taxonomic profiling completed."