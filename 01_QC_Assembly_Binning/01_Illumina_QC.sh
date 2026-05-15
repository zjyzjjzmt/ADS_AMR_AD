#!/bin/bash
# Description: Quality control for Illumina paired-end reads and sequence statistics.

# Activate environment
source deactivate
conda activate metawrap-env

# Create necessary directories
mkdir -p 01READ_QC 01Reads_fq

echo "[INFO] Starting MetaWRAP read_qc..."
for F in 01RAW_READS/*_1.fastq.gz; do 
    R="${F%_1.fastq.gz}_2.fastq.gz"
    BASE=$(basename "$F")
    SAMPLE="${BASE%_1.fastq.gz}"
    
    if [ ! -d "01READ_QC/$SAMPLE" ]; then
        # Parameters: -t (threads), --skip-bmtagger (skip human sequence removal if not applicable)
        metawrap read_qc -1 "$F" -2 "$R" -t 20 --skip-bmtagger -o "01READ_QC/$SAMPLE"
    else
        echo "[WARN] 01READ_QC/$SAMPLE already exists, skipping QC..."
    fi
done

echo "[INFO] Moving clean reads to centralized folder..."
for i in 01READ_QC/*; do 
    if [ -d "$i" ]; then  
        SAMPLE=$(basename "$i")
        R1_FILE="${i}/${SAMPLE}_1_val_1.fq.gz"
        R2_FILE="${i}/${SAMPLE}_2_val_2.fq.gz"
        
        if [ -f "$R1_FILE" ] && [ -f "$R2_FILE" ]; then
            mv "$R1_FILE" "01Reads_fq/${SAMPLE}_1.fastq.gz"
            mv "$R2_FILE" "01Reads_fq/${SAMPLE}_2.fastq.gz"
            echo "[INFO] $SAMPLE transfer completed."
        fi
    fi
done

echo "[INFO] Generating read counts statistics..."
STATS_FILE="01Reads_fq/Clean_reads_number.txt"
echo "Sample_Name: Reads_Count" > "$STATS_FILE"

for file in 01Reads_fq/*.fastq.gz; do
    if [ -f "$file" ]; then
        lines=$(gzip -dc "$file" | wc -l)
        reads=$((lines / 4))
        echo "$(basename "$file" .fastq.gz): ${reads} reads" >> "$STATS_FILE"
    fi
done
echo "[INFO] Pipeline finished successfully."