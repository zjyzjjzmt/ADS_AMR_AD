#!/bin/bash


source deactivate
conda activate metawrap-env

for F in 01Reads_fq/*_1.fastq.gz; do
    R=${F%_*}_2.fastq.gz
    BASE=${F##*/}
    SAMPLE=${BASE%_1.fastq.gz} 
    
    echo "========== Processing Sample: $SAMPLE =========="

    # 1. Assembly using MEGAHIT via metaWRAP
    # Parameter -m 500 sets minimum contig length to 500bp
    echo "[INFO] Starting Assembly..."
    metawrap assembly -1 $F -2 $R -m 500 -t 50 -o 05ASSEMBLY/$SAMPLE
    rm -rf 05ASSEMBLY/$SAMPLE/megahit

    # 2. Temporary decompression for Binning module requirements
    echo "[INFO] Temporarily extracting reads..."
    TMP_R1="01Reads_fq/tmp_${SAMPLE}_1.fastq"
    TMP_R2="01Reads_fq/tmp_${SAMPLE}_2.fastq"
    zcat $F > $TMP_R1
    zcat $R > $TMP_R2

    # 3. Integrated Binning (MetaBAT2, CONCOCT, MaxBin2)
    echo "[INFO] Starting Binning..."
    metawrap binning -o 06MAGs/01bins/$SAMPLE -t 50 -m 500 \
        -a 05ASSEMBLY/$SAMPLE/final_assembly.fasta \
        --metabat2 --concoct --maxbin2 $TMP_R1 $TMP_R2 --universal
    rm -rf 06MAGs/01bins/$SAMPLE/work_files

    # 4. Bin Refinement
    # Parameters: -c 50 (Minimum 50% completeness) -x 20 (Maximum 20% contamination)
    # Note: Further stringent filtering (>80% comp, <5% cont) will be applied in downstream analysis.
    echo "[INFO] Starting Bin Refinement..."
    metawrap bin_refinement -o 06MAGs/02refine_bins/$SAMPLE -t 50 \
        -A 06MAGs/01bins/$SAMPLE/maxbin2_bins/ \
        -B 06MAGs/01bins/$SAMPLE/metabat2_bins/ \
        -C 06MAGs/01bins/$SAMPLE/concoct_bins/ \
        -c 50 -x 20
    rm -rf 06MAGs/02refine_bins/$SAMPLE/work_files

    # 5. Clean up temporary files
    echo "[INFO] Cleaning temporary files..."
    rm $TMP_R1 $TMP_R2

    echo "========== Sample $SAMPLE Finished =========="
done