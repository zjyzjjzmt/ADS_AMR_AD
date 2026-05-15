#!/bin/bash
#SBATCH --job-name=viral_taxonomy
#SBATCH --cpus-per-task=50

source deactivate

# ==========================================
# 1. Classification via geNomad
# ==========================================
echo "[INFO] Running geNomad taxonomy..."
conda activate genomad
mkdir -p 02genomad
genomad end-to-end --cleanup 01Clustergenome/DWTP_VCs.fasta 02genomad /path/to/genomad_db

# ==========================================
# 2. Classification via CAT (Contig Annotation Tool)
# ==========================================
echo "[INFO] Running CAT taxonomy..."
conda deactivate
conda activate tools
mkdir -p 02CAT

CAT contigs -c 01Clustergenome/DWTP_VCs.fasta \
            -d /path/to/CAT_database \
            -t /path/to/CAT_taxonomy
            
CAT add_names -i out.CAT.contig2classification.txt -o CAT_contigs_classification.txt -t /path/to/CAT_taxonomy
CAT add_names -i out.CAT.ORF2LCA.txt -o CAT_ORFs_classification.txt -t /path/to/CAT_taxonomy

mv out.CAT.* 02CAT/
mv CAT_contigs_classification.txt 02CAT/
mv CAT_ORFs_classification.txt 02CAT/
rm -rf out.CAT.alignment.diamond

# ==========================================
# 3. Network-based Classification via vConTACT2 (For contigs > 10kb)
# ==========================================
echo "[INFO] Running vConTACT2..."
conda deactivate
conda activate vcontact2
mkdir -p 04vcontact2

seqtk seq -L 10000 01Clustergenome/DWTP_VCs.fasta > 04vcontact2/DWTP_VCs_10kb.fasta
prodigal -i 04vcontact2/DWTP_VCs_10kb.fasta -o 04vcontact2/DWTP_VCs_10kb.genes -a 04vcontact2/DWTP_VCs_10kb.faa -p meta

vcontact2_gene2genome -p 04vcontact2/DWTP_VCs_10kb.faa -o 04vcontact2/DWTP_VCs_10kb_g2g.csv -s 'Prodigal-FAA'
vcontact2 -r 04vcontact2/DWTP_VCs_10kb.faa \
          --rel-mode 'Diamond' \
          --proteins-fp 04vcontact2/DWTP_VCs_10kb_g2g.csv \
          --db 'ProkaryoticViralRefSeq201-Merged' \
          --pcs-mode MCL \
          --vcs-mode ClusterONE \
          --c1-bin /path/to/cluster_one-1.0.jar \
          --output-dir 04vcontact2/vContactOut/

# ==========================================
# 4. Custom Majority-Rules Classification against Viral RefSeq
# ==========================================
echo "[INFO] Running Majority-Rules custom pipeline..."
conda deactivate
conda activate tools
mkdir -p 08majority_rules/02blastp_accession

# Protein alignment against Viral RefSeq
diamond blastp -q 01Clustergenome/DWTP_VCs.faa -d /path/to/viral_refseq/viral.1.protein.faa \
               -o 08majority_rules/DWTP_VCs_blastp.txt \
               --query-cover 50 --subject-cover 50 --evalue 1e-5 -k 1

# Filter hits by bitscore (>=50) and extract accession IDs
csvtk filter -t -f "12>=50" 08majority_rules/DWTP_VCs_blastp.txt > 08majority_rules/DWTP_VCs_blastp_50score.txt
cut -f 2 08majority_rules/DWTP_VCs_blastp_50score.txt > 08majority_rules/DWTP_VCs_blastp_accession.txt

# Split for parallel mapping to Taxonomy IDs
split -l 1000 08majority_rules/DWTP_VCs_blastp_accession.txt 08majority_rules/02blastp_accession/part_
for file in 08majority_rules/02blastp_accession/part_*; do
    rg -f "$file" /path/to/prot.accession2taxid.FULL --no-line-number >> 08majority_rules/DWTP_VCs_refseq_accession2taxid.txt
done

# Extract lineage using TaxonKit and aggregate taxonomy to Family level using Awk
cut -f 2 08majority_rules/DWTP_VCs_refseq_accession2taxid.txt > 08majority_rules/DWTP_VCs_blastp_taxid.txt
taxonkit lineage 08majority_rules/DWTP_VCs_blastp_taxid.txt --data-dir /path/to/taxdump | taxonkit reformat -F --data-dir /path/to/taxdump | cut -f 1,3 > 08majority_rules/DWTP_VCs_taxid_taxonomy.txt
paste 08majority_rules/DWTP_VCs_refseq_accession2taxid.txt 08majority_rules/DWTP_VCs_taxid_taxonomy.txt | cut -f 1,4 > 08majority_rules/DWTP_VCs_accession_taxonomy.txt

# Map taxonomy back to original hits and format string
awk -F'\t' 'FNR==NR{a[$1]=$2; next}; {if($2 in a) {print $0, "\t"a[$2];} else {print $0, "\tNA"}}' 08majority_rules/DWTP_VCs_accession_taxonomy.txt 08majority_rules/DWTP_VCs_blastp_50score.txt > 08majority_rules/DWTP_VCs_blastp_viref_50_tax.txt
awk -F'\t' -v OFS='\t' '{sub(/_[0-9]+$/, "", $1); split($13, a, ";"); $13=a[1] ";" a[2] ";" a[3] ";" a[4] ";" a[5]; print $1, $13, 1}' 08majority_rules/DWTP_VCs_blastp_viref_50_tax.txt > 08majority_rules/DWTP_VCs_family.txt

# Aggregate counts to define the majority taxonomy per contig
awk -F'\t' '{sums[$1"\t"$2] += $3} END {for (key in sums) print key"\t"sums[key]}' 08majority_rules/DWTP_VCs_family.txt > 08majority_rules/DWTP_VCs_family_stax.txt
cut -f 1,3 08majority_rules/DWTP_VCs_family_stax.txt | datamash -sW -g1 sum 2 > 08majority_rules/DWTP_VCs_contigs_stax.txt
awk 'FNR==NR{a[$1]=$2; next}; {if($1 in a) {print $0, "\t"a[$1];} else {print $0, "\tNA"}}' 08majority_rules/DWTP_VCs_contigs_stax.txt 08majority_rules/DWTP_VCs_family_stax.txt > 08majority_rules/DWTP_VCs_family_contigs.txt

# ==========================================
# 5. Database Matching & Host Information (IMG/VR & RefSeq)
# ==========================================
echo "[INFO] Searching against known viral databases..."
mkdir -p 09blast2IMGVR

# Blast against RefSeq Genomic
blastn -query 01Clustergenome/DWTP_VCs.fasta -db /path/to/viral.1.1.genomic.fna -out 09blast2IMGVR/DWTP_VCs_blastnVG.txt \
       -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident staxids sscinames scomnames sblastnames bitscore salltitles qcovs qcovhsp stitle" \
       -dust no -max_target_seqs 1 -perc_identity 95 -evalue 0.00001 -num_threads 20

# Blast against IMG/VR and extract host context
blastn -query 01Clustergenome/DWTP_VCs.fasta -db /path/to/IMGVR -out 09blast2IMGVR/DWTP_VCs_blastnIMGVR.txt \
       -outfmt "6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident staxids sscinames scomnames sblastnames bitscore salltitles qcovs qcovhsp stitle" \
       -dust no -max_target_seqs 1 -perc_identity 95 -evalue 0.00001 -num_threads 20
       
cut -f 2 09blast2IMGVR/DWTP_VCs_blastnIMGVR.txt > 09blast2IMGVR/DWTP_IMGVR_VCblastn.txt
cat /path/to/IMGVR_all_Host_information-high_confidence.tsv | csvtk grep -t -P 09blast2IMGVR/DWTP_IMGVR_VCblastn.txt > 09blast2IMGVR/IMGVR_host_results.txt

# ==========================================
# 6. Novelty Assessment at the Protein Level
# ==========================================
echo "[INFO] Assessing viral novelty..."
mkdir -p 10novelty

# Parameters: id 30%, cov 50%, evalue 1e-5
diamond blastp -q 01Clustergenome/DWTP_VCs.faa -d /path/to/IMGVR_proteins -o 10novelty/01DWTP_VCs_IMGVR.txt --id 30 --query-cover 50 --evalue 1e-5 -k 1
diamond blastp -q 01Clustergenome/DWTP_VCs.faa -d /path/to/viral_refseq -o 10novelty/02DWTP_VCs_Refseq.txt --id 30 --query-cover 50 --evalue 1e-5 -k 1

echo "[INFO] Taxonomic and novelty assessments completed."