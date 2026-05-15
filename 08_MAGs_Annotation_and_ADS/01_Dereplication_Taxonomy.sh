#!/bin/bash
#SBATCH --job-name=mag_drep_gtdbtk
#SBATCH --cpus-per-task=50

source deactivate

# 1. Prepare Bins for dRep
echo "[INFO] Aggregating and renaming bins..."
conda activate tools
mkdir -p 06MAGs/03beforedRep 06MAGs/04dRep

for SAMPLE_DIR in 06MAGs/02refine_bins/*; do
    SAMPLE=$(basename "$SAMPLE_DIR")
    BIN_DIR="$SAMPLE_DIR/metawrap_50_20_bins"
    
    if [ -d "$BIN_DIR" ]; then
        for fa in "$BIN_DIR"/*.fa; do
            # Add sample prefix if not already present
            if [[ $(basename "$fa") != $SAMPLE-* ]]; then
                mv "$fa" "$BIN_DIR/${SAMPLE}-$(basename "$fa")"
            fi
        done
        cp "$BIN_DIR"/*.fa 06MAGs/03beforedRep/
    fi
done

# 2. De-replicate Genomes (Completeness >70%, Contamination <10%)
echo "[INFO] Running dRep..."
conda deactivate
conda activate dRep

dRep dereplicate 06MAGs/04dRep/ -g 06MAGs/03beforedRep/*.fa -comp 70 -con 10 -p 50

# 3. Taxonomic Classification using GTDB-Tk
echo "[INFO] Running GTDB-Tk..."
conda deactivate
conda activate gtdbtk
mkdir -p 06MAGs/06gtdbtk_classify

gtdbtk classify_wf --skip_ani_screen \
                   --genome_dir 06MAGs/04dRep/dereplicated_genomes/ \
                   -x fa \
                   --out_dir 06MAGs/06gtdbtk_classify/ \
                   --cpus 50