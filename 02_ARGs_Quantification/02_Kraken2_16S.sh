#!/bin/bash

source deactivate
conda activate tools

# Path to Kraken2 16S database
K2_DB="/path/to/kraken2/16S_Greengenes_k2db"
mkdir -p 04reads_based_analysis/07kraken2

echo "[INFO] Running Kraken2 for 16S extraction..."
for i in 01Reads_fq/*.fastq.gz; do
    kraken2 --db $K2_DB "$i" \
            --output "$i.txt" \
            --report "$i.report.txt" \
            --classified-out "$i.16s.fasta" \
            --threads 60
            
    mv "$i.report.txt" 04reads_based_analysis/07kraken2/
    mv "$i.16s.fasta" 04reads_based_analysis/07kraken2/
    mv "$i.txt" 04reads_based_analysis/07kraken2/
done

# Count the number of 16S reads for downstream normalization
cd 04reads_based_analysis/07kraken2
echo "Sample_Name: 16S_Reads_Count" > 16s_reads_number.txt
for file in *.fasta; do
    # Note: Using ^> to accurately match FASTA sequence headers
    echo "$(basename "$file" .fasta): $(grep -c "^>" "$file") reads" >> 16s_reads_number.txt
done
echo "[INFO] 16S read extraction finished."