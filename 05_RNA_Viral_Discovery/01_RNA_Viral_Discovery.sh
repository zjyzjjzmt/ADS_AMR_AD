#!/bin/bash
#SBATCH --job-name=rna_viral_discovery
#SBATCH --cpus-per-task=50

source deactivate

mkdir -p 08RNA_VCs/00Contigs 08RNA_VCs/rdrp/test 08RNA_VCs/01genomad

# 1. ORF Prediction on >1000bp contigs
conda activate tools
prodigal -i 08RNA_VCs/00Contigs/AD_RNA_contigs1000bp.fasta \
         -o 08RNA_VCs/00Contigs/AD_RNA_contigs1000bp.genes \
         -a 08RNA_VCs/00Contigs/AD_RNA_contigs1000bp.faa \
         -d 08RNA_VCs/00Contigs/AD_RNA_contigs1000bp.fna -p meta
cp 08RNA_VCs/00Contigs/AD_RNA_contigs1000bp.faa 08RNA_VCs/rdrp/test/

# 2. RdRp-based Identification (HMMER)
echo "[INFO] Running RdRp-search..."
conda deactivate
conda activate rdrp
cd 08RNA_VCs/rdrp
./rdrp-search -f "test/*.faa" -x 4 -s 50 -t 1 -g -o test/out
cd ../../

# 3. geNomad Identification
echo "[INFO] Running geNomad..."
conda deactivate
conda activate genomad
genomad end-to-end --cleanup 08RNA_VCs/00Contigs/AD_RNA_contigs1000bp.fasta \
                   08RNA_VCs/01genomad /path/to/genomad_db

# 4. Merge Results (Riboviria + RdRp hits) and De-duplicate
echo "[INFO] Merging candidates..."
conda deactivate
conda activate tools

# Extract Riboviria from geNomad
grep "Riboviria" 08RNA_VCs/01genomad/AD_RNA_contigs1000bp_summary/AD_RNA_contigs1000bp_virus_summary.tsv | cut -f 1 > 08RNA_VCs/RNA_virus.txt
# Extract hits from RdRp search and clean sequence IDs
awk '{sub(/^[^|]*\|\|/, "", $1); sub(/_[^_]*$/, "", $1); print $1}' 08RNA_VCs/rdrp/test/out/final-hits.fa.fai >> 08RNA_VCs/RNA_virus.txt

sort 08RNA_VCs/RNA_virus.txt | uniq > 08RNA_VCs/RNA_virus_nr.txt
seqtk subseq 08RNA_VCs/00Contigs/AD_RNA_contigs1000bp.fasta 08RNA_VCs/RNA_virus_nr.txt > 08RNA_VCs/RNA_virus.fasta

# CD-HIT De-duplication (100% identity)
cd-hit-est -i 08RNA_VCs/RNA_virus.fasta -o 08RNA_VCs/AD_RNA_virus_nr.fasta -c 1.0 -n 10 -M 16000 -T 8