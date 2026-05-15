#!/bin/bash
#SBATCH --job-name=hic_mapping
#SBATCH --cpus-per-task=50

source deactivate
conda activate tools

READS_DIR="01Reads_fq"
OUTPUT_DIR="02Hi-C"
# Reference must contain both MAGs (bins) and viral contigs
ASSEMBLY="hic_assembly.fasta"

mkdir -p "$OUTPUT_DIR"

# 1. Build BWA index
echo "[INFO] Building BWA index for combined reference..."
bwa index "$ASSEMBLY"

SAMPLE_COUNT=1

for R1 in "$READS_DIR"/*_1.fastq; do
    SAMPLE=$(basename "$R1" _1.fastq)
    R2="${READS_DIR}/${SAMPLE}_2.fastq"

    if [[ -f "$R2" ]]; then
        echo "========== Processing sample: $SAMPLE (ID: $SAMPLE_COUNT) =========="

        # 2. BWA MEM Mapping (Critical: -5SP for Hi-C data)
        bwa mem -5SP -t 50 "$ASSEMBLY" "$R1" "$R2" > "$OUTPUT_DIR/${SAMPLE}.sam"

        # 3. Filter and Sort BAM (Critical: -F 0x904)
        samtools view -b -F 0x904 -@ 50 "$OUTPUT_DIR/${SAMPLE}.sam" > "$OUTPUT_DIR/${SAMPLE}_FILTERED.bam"
        samtools sort -n -@ 50 "$OUTPUT_DIR/${SAMPLE}_FILTERED.bam" -o "$OUTPUT_DIR/${SAMPLE}_NAME_SORTED.bam"
        rm "$OUTPUT_DIR/${SAMPLE}.sam"

        # 4. Convert BAM to BEDPE format to extract cross-links
        bedtools bamtobed -bedpe -i "$OUTPUT_DIR/${SAMPLE}_NAME_SORTED.bam" > "$OUTPUT_DIR/${SAMPLE}.bedpe"

        # 5. Extract strictly cross-contig links (exclude intra-contig)
        awk '$1 != $4' "$OUTPUT_DIR/${SAMPLE}.bedpe" > "$OUTPUT_DIR/${SAMPLE}_FILTERED.bedpe"

        # 6. Count interaction frequencies between contig pairs
        awk '{print $1 "\t" $4}' "$OUTPUT_DIR/${SAMPLE}_FILTERED.bedpe" | sort | uniq -c > "$OUTPUT_DIR/${SAMPLE}_contig_links.txt"

        # 7. Format output (Count \t Contig_A \t Contig_B)
        awk '{print $1 "\t" $2 "\t" $3}' "$OUTPUT_DIR/${SAMPLE}_contig_links.txt" > "$OUTPUT_DIR/${SAMPLE}_contig_links_modified.txt"

        # 8. Strict filtering: Keep only interactions with > 5 links
        awk '$1 > 5' "$OUTPUT_DIR/${SAMPLE}_contig_links_modified.txt" > "$OUTPUT_DIR/${SAMPLE}_contig_links_step1.txt"

        # 9. Exclude internal MAG-MAG (bin-bin) linkages
        awk '!($2 ~ /bin/ && $3 ~ /bin/)' "$OUTPUT_DIR/${SAMPLE}_contig_links_step1.txt" > "$OUTPUT_DIR/${SAMPLE}_contig_links_step2.txt"

        # 10. Isolate Phage-Host interactions (Must involve exactly one bin)
        awk '($2 ~ /bin/ || $3 ~ /bin/)' "$OUTPUT_DIR/${SAMPLE}_contig_links_step2.txt" > "$OUTPUT_DIR/${SAMPLE}_contig_links_final.txt"

        # 11. Clean bin names (Remove trailing identifiers added by binners)
        awk '{
            if ($2 ~ /bin/) { gsub(/_.*$/, "", $2) }
            if ($3 ~ /bin/) { gsub(/_.*$/, "", $3) }
            print $1 "\t" $2 "\t" $3
        }' "$OUTPUT_DIR/${SAMPLE}_contig_links_final.txt" > "$OUTPUT_DIR/${SAMPLE}_contig_links_final_modified.txt"

        # 12. Aggregate and deduplicate identical pairs
        awk '{
            key = $2 "\t" $3
            sum[key] += $1
        }
        END {
            for (key in sum) {
                print sum[key] "\t" key
            }
        }' "$OUTPUT_DIR/${SAMPLE}_contig_links_final_modified.txt" > "$OUTPUT_DIR/${SAMPLE}_contig_links_no_duplicates.txt"

        echo "[INFO] Finished sample: $SAMPLE"
        SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    fi
done

echo "[INFO] Hi-C Phage-Host mapping completely finished!"