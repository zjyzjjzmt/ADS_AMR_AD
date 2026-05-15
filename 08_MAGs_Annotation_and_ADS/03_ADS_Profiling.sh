#!/bin/bash
#SBATCH --job-name=ads_profiling
#SBATCH --cpus-per-task=40

source deactivate
MAGS_DIR="06MAGs/04dRep/dereplicated_genomes"
PROKKA_DIR="06MAGs/07Prokka"

# 1. Structural Annotation via Prokka (Required format for PADLOC/DefenseFinder)
echo "[INFO] Running Prokka for structural annotation..."
conda activate tools
for F in ${MAGS_DIR}/*.fa; do
    SAMPLE=$(basename "$F" .fa)
    prokka --outdir $PROKKA_DIR/$SAMPLE --prefix $SAMPLE --cpus 40 --metagenome $F
done

# 2. ADS Mining using PADLOC
echo "[INFO] Running PADLOC..."
conda deactivate
conda activate padloc
mkdir -p 06MAGs/Defense_Systems/01PADLOC

for DIR in $PROKKA_DIR/*; do
    SAMPLE=$(basename "$DIR")
    padloc --faa $DIR/$SAMPLE.faa \
           --gff $DIR/$SAMPLE.gff \
           --outdir 06MAGs/Defense_Systems/01PADLOC/$SAMPLE \
           --cpu 40
done

# 3. ADS Mining using DefenseFinder
echo "[INFO] Running DefenseFinder..."
conda deactivate
conda activate defensefinder
mkdir -p 06MAGs/Defense_Systems/02DefenseFinder

for DIR in $PROKKA_DIR/*; do
    SAMPLE=$(basename "$DIR")
    defense-finder run $DIR/$SAMPLE.faa \
                   --out-dir 06MAGs/Defense_Systems/02DefenseFinder/$SAMPLE
done

echo "[INFO] ADS profiling and entire analytical pipeline successfully completed."