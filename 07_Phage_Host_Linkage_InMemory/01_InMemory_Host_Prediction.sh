#!/bin/bash
#SBATCH --job-name=in_memory_host_linkage
#SBATCH --cpus-per-task=40

source deactivate

# Define central paths
MAGS_DIR="06MAGs/04dRep/dereplicated_genomes"
VCS_FASTA="01Clustergenome/DWTP_VCs.fasta"
OUT_DIR="10host-link"

mkdir -p ${OUT_DIR}/01crispr ${OUT_DIR}/02tRNA ${OUT_DIR}/03homology 06MAGs/04dRep/01MAGs_crisprcas 06MAGs/04dRep/03MAGs_spacers

# ==========================================
# 1. CRISPR-Cas Spacer Extraction from MAGs
# ==========================================
echo "[INFO] Extracting CRISPR spacers using CRT..."
conda activate tools

for F in ${MAGS_DIR}/*.fa; do
    BASE=${F##*/}
    SAMPLE=${BASE%.fa*}
    
    # Run CRISPR Recognition Tool (CRT)
    java -cp CRT1.2-CLI.jar crt $F $F.out
    
    # Standardize contig headers to include sample/MAG ID
    sed -i "s/k141\|k119/$SAMPLE/g" $F.out
    mv $F.out 06MAGs/04dRep/01MAGs_crisprcas/
done

# Extract spacers using custom R script and format FASTA
cd 06MAGs/04dRep
Rscript MAGs_CRT_spacers.R
cat 03MAGs_spacers/*.txt > MAGs_spacers.fasta
seqtk seq -L 1 MAGs_spacers.fasta > ../../${OUT_DIR}/01crispr/MAGs_spacers_nr.fasta
cd ../../

# ==========================================
# 2. tRNA Extraction from Viral Contigs
# ==========================================
echo "[INFO] Extracting viral tRNAs using Aragorn..."
aragorn -t -fasta -wa -o ${OUT_DIR}/02tRNA/DWTP_VCs_tRNA.fasta ${VCS_FASTA}

# Parse and clean tRNA outputs
python3 ${OUT_DIR}/02tRNA/tRNA.py

# ==========================================
# 3. Centralized MAG Database Preparation
# ==========================================
echo "[INFO] Building centralized MAG database..."
mkdir -p 06MAGs/04dRep/02MAGs_database

for F in ${MAGS_DIR}/*.fa; do
    BASE=${F##*/}
    SAMPLE=${BASE%.fa*}
    sed "s/k141\|k119/$SAMPLE/g" $F > 06MAGs/04dRep/02MAGs_database/$SAMPLE.fa
done

cat 06MAGs/04dRep/02MAGs_database/*.fa > ${OUT_DIR}/02tRNA/DWTP_MAGs_database.fa
cp ${OUT_DIR}/02tRNA/DWTP_MAGs_database.fa ${OUT_DIR}/03homology/

# Build BLAST databases
makeblastdb -dbtype nucl -in ${OUT_DIR}/02tRNA/DWTP_MAGs_database.fa -input_type fasta -title DWTP_MAGs -out ${OUT_DIR}/02tRNA/DWTP_MAGs -blastdb_version 5
makeblastdb -dbtype nucl -in ${VCS_FASTA} -input_type fasta -title DWTP_VCs -out ${OUT_DIR}/01crispr/DWTP_VCs -blastdb_version 5

# ==========================================
# 4. Tri-Evidence Linkage (BLAST Mapping)
# ==========================================
echo "[INFO] Running Tri-Evidence Mapping..."

# A. tRNA Linkage (Viral tRNA against MAGs | Strict 100% ID)
echo " -> Mapping tRNAs..."
blastn -query ${OUT_DIR}/02tRNA/DWTP_VCs_tRNA_nr.fasta \
       -db ${OUT_DIR}/02tRNA/DWTP_MAGs \
       -out ${OUT_DIR}/02tRNA/DWTP_MAGs_hosts_tRNA.txt \
       -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident staxids bitscore salltitles qcovs qcovhsp stitle" \
       -dust no -perc_identity 100 -evalue 0.0001 -num_threads 40 

# B. CRISPR Linkage (MAG Spacers against Viral Contigs | 97% ID, short-task)
echo " -> Mapping CRISPR Spacers..."
blastn -task blastn-short \
       -query ${OUT_DIR}/01crispr/MAGs_spacers_nr.fasta \
       -db ${OUT_DIR}/01crispr/DWTP_VCs \
       -out ${OUT_DIR}/01crispr/DWTP_VCs_MAGs_criprcas_host.txt \
       -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident bitscore mismatch qcovs qcovhsp" \
       -perc_identity 97 -max_target_seqs 1 -num_threads 40

# C. Homology Linkage (Whole Viral Contigs against MAGs | 80% ID)
echo " -> Mapping Sequence Homology..."
blastn -query ${VCS_FASTA} \
       -db ${OUT_DIR}/03homology/DWTP_MAGs \
       -out ${OUT_DIR}/03homology/DWTP_Virus_homology.txt \
       -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident staxids sscinames scomnames sblastnames bitscore salltitles qcovs qcovhsp stitle" \
       -dust no -perc_identity 80 -evalue 0.00001 -num_threads 40

echo "[INFO] In-Memory Phage-Host Linkage Analysis Completed."