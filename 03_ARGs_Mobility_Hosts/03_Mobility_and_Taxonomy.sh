#!/bin/bash
#SBATCH --job-name=ARC_mobility_taxonomy

cd 07ARCs_analysis

# 1. Plasmid Prediction using PlasFlow
echo "[INFO] Running PlasFlow..."
source deactivate
conda activate plasflow

mkdir -p 08plasflow
PlasFlow.py --input DWTP_ARCs_nr.fasta \
            --output 08plasflow/DWTP_ARCs_nr.plasflow.txt \
            --threshold 0.7

# 2. Taxonomic Classification using CAT (Contig Annotation Tool)
echo "[INFO] Running CAT classification..."
source deactivate
conda activate tools
mkdir -p 07CAT

CAT contigs -c DWTP_ARCs_nr.fasta \
            -d /path/to/CAT_database \
            -t /path/to/CAT_taxonomy \
            -p DWTP_ARCs_nr.pro.fa
            
CAT add_names -i out.CAT.contig2classification.txt -o 07CAT/CAT_contigs_classification.txt -t /path/to/CAT_taxonomy
CAT add_names -i out.CAT.ORF2LCA.txt -o 07CAT/CAT_ORFs_classification.txt -t /path/to/CAT_taxonomy

mv out.CAT.* 07CAT/
rm -rf out.CAT.alignment.diamond

# 3. Taxonomic Classification using Taxator-tk (via MetaWRAP)
echo "[INFO] Running Taxator-tk via MetaWRAP..."
source deactivate
conda activate metawrap-env
mkdir -p 10contigs
cp DWTP_ARCs_nr.fasta 10contigs/

metawrap classify_bins -b 10contigs -o 11metawrap-classify -t 10

echo "[INFO] Mobility and taxonomy analyses successfully completed."