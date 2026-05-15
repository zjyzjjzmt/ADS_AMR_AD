#!/bin/bash
#SBATCH --job-name=mag_functional_annotation

source deactivate
conda activate tools

MAGS_DIR="06MAGs/04dRep/dereplicated_genomes"
mkdir -p 06MAGs/ORFs_genes 06MAGs/ORFs_pro 06MAGs/ORFs_nucl 06MAGs/06RNA

# 1. ORF Prediction & rRNA extraction
echo "[INFO] Predicting ORFs and rRNAs..."
for F in ${MAGS_DIR}/*.fa; do
    SAMPLE=$(basename "$F" .fa)
    
    # Prodigal
    prodigal -i $F -o 06MAGs/ORFs_genes/$SAMPLE.genes \
             -a 06MAGs/ORFs_pro/$SAMPLE.pro.fa \
             -d 06MAGs/ORFs_nucl/$SAMPLE.fa -p meta
             
    # Barrnap (rRNA)
    barrnap --outseq 06MAGs/06RNA/$SAMPLE.RNA.fa $F
done

# 2. Prepare Databases
echo "[INFO] Formatting DIAMOND and BLAST databases..."
DB_DIR="path/to/databases"
diamond makedb --in $DB_DIR/SARG.2.2.fasta -d 06MAGs/ORFs_nucl/01SARG/SARG.2.2_nr
diamond makedb --in $DB_DIR/CARD3.1.4.fasta -d 06MAGs/ORFs_nucl/02card/card3.1.4_nr
diamond makedb --in $DB_DIR/BacMet2.fasta -d 06MAGs/ORFs_nucl/03bacmet/bacmet2_nr
diamond makedb --in $DB_DIR/ICE3.0.fasta -d 06MAGs/ORFs_nucl/04ICEs/ICE3.0_nr
diamond makedb --in $DB_DIR/victors_pro.fasta -d 06MAGs/ORFs_nucl/06victors/victors_nr
makeblastdb -dbtype nucl -in $DB_DIR/MGEs.fasta -title MGEs -out 06MAGs/ORFs_nucl/05MGEs/MGEs 

# 3. Functional Mapping
echo "[INFO] Running functional alignments..."
for F in 06MAGs/ORFs_nucl/*.fa; do
    # DIAMOND blastx (Strict parameters)
    diamond blastx -d 06MAGs/ORFs_nucl/01SARG/SARG.2.2_nr -q $F -o $F-SARG.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
    diamond blastx -d 06MAGs/ORFs_nucl/02card/card3.1.4_nr -q $F -o $F-card.txt --evalue 1e-5 --query-cover 70 --id 80 -k 1
    diamond blastx -d 06MAGs/ORFs_nucl/03bacmet/bacmet2_nr -q $F -o $F-bacmet2.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
    diamond blastx -d 06MAGs/ORFs_nucl/04ICEs/ICE3.0_nr -q $F -o $F-ICE.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
    diamond blastx -d 06MAGs/ORFs_nucl/06victors/victors_nr -q $F -o $F-victors.txt --evalue 1e-10 --query-cover 70 --id 80 -k 1
    
    # BLASTN for MGEs
    blastn -query $F -db 06MAGs/ORFs_nucl/05MGEs/MGEs -out $F-MGEs.txt -evalue 1e-10 -perc_identity 80 -num_threads 10 -outfmt 6 -max_target_seqs 1
    
    # Organize outputs
    mv $F-SARG.txt 06MAGs/ORFs_nucl/01SARG/
    mv $F-card.txt 06MAGs/ORFs_nucl/02card/
    mv $F-bacmet2.txt 06MAGs/ORFs_nucl/03bacmet/
    mv $F-ICE.txt 06MAGs/ORFs_nucl/04ICEs/
    mv $F-victors.txt 06MAGs/ORFs_nucl/06victors/
    mv $F-MGEs.txt 06MAGs/ORFs_nucl/05MGEs/
done