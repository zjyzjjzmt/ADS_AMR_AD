#!/bin/bash
#SBATCH --job-name=ARC_extraction
#SBATCH --time=1-00:00:00

source deactivate
conda activate tools

mkdir -p 05ASSEMBLY/ORFs_genes 05ASSEMBLY/ORFs_protein 05ASSEMBLY/ORFs_nucl
mkdir -p 05ASSEMBLY/ORFs_nucl/01SARG 05ASSEMBLY/ORFs_nucl/02card

echo "[INFO] Predicting ORFs using Prodigal..."
for asm in 05ASSEMBLY/*/final_assembly.fasta; do
    SAMPLE=$(basename $(dirname "$asm"))
    prodigal -i "$asm" \
             -o 05ASSEMBLY/$SAMPLE/$SAMPLE.genes \
             -a 05ASSEMBLY/$SAMPLE/$SAMPLE.pro.fa \
             -d 05ASSEMBLY/$SAMPLE/$SAMPLE.fa -p meta
             
    mv 05ASSEMBLY/$SAMPLE/$SAMPLE.genes 05ASSEMBLY/ORFs_genes/
    mv 05ASSEMBLY/$SAMPLE/$SAMPLE.pro.fa 05ASSEMBLY/ORFs_protein/
    mv 05ASSEMBLY/$SAMPLE/$SAMPLE.fa 05ASSEMBLY/ORFs_nucl/
done

echo "[INFO] Extracting ARCs based on SARG & CARD databases..."
# Paths to databases
SARG_DB="path/to/db/SARG.2.2_nr"
CARD_DB="path/to/db/card3.1.4_nr"

for fa in 05ASSEMBLY/ORFs_nucl/*.fa; do
    SAMPLE=$(basename "$fa" .fa)
    
    # 1. SARG Annotation
    diamond blastx -d $SARG_DB -q "$fa" -o "05ASSEMBLY/ORFs_nucl/01SARG/${SAMPLE}.fa.txt" --evalue 1e-10 --query-cover 70 --id 80 -k 1
    # Parse Prodigal headers to extract source contig IDs
    cut -f 1,2 "05ASSEMBLY/ORFs_nucl/01SARG/${SAMPLE}.fa.txt" | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' | cut -f 1,2,3,4,5,6,8 | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' | sort -k 1n > "05ASSEMBLY/ORFs_nucl/01SARG/${SAMPLE}_ARGs.txt"
    cut -d' ' -f1 "05ASSEMBLY/ORFs_nucl/01SARG/${SAMPLE}_ARGs.txt" | sort -u > "05ASSEMBLY/ORFs_nucl/01SARG/${SAMPLE}_ARGs_list.txt"
    # Extract full contig sequences using seqtk
    seqtk subseq "05ASSEMBLY/$SAMPLE/final_assembly.fasta" "05ASSEMBLY/ORFs_nucl/01SARG/${SAMPLE}_ARGs_list.txt" | sed "s/k141/$SAMPLE/g" > "05ASSEMBLY/ORFs_nucl/01SARG/${SAMPLE}_ARGs.fa"

    # 2. CARD Annotation
    diamond blastx -d $CARD_DB -q "$fa" -o "05ASSEMBLY/ORFs_nucl/02card/${SAMPLE}.fa.txt" --evalue 1e-5 --query-cover 70 --id 80 -k 1
    # Parse Prodigal headers
    cut -f 1,2 "05ASSEMBLY/ORFs_nucl/02card/${SAMPLE}.fa.txt" | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' | cut -f 1,2,3,4,5,6,8 | awk -F " " '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6" "$7}' | sort -k 1n > "05ASSEMBLY/ORFs_nucl/02card/${SAMPLE}_ARGs.txt"
    cut -d' ' -f1 "05ASSEMBLY/ORFs_nucl/02card/${SAMPLE}_ARGs.txt" | sort -u > "05ASSEMBLY/ORFs_nucl/02card/${SAMPLE}_ARGs_list.txt"
    seqtk subseq "05ASSEMBLY/$SAMPLE/final_assembly.fasta" "05ASSEMBLY/ORFs_nucl/02card/${SAMPLE}_ARGs_list.txt" | sed "s/k141/$SAMPLE/g" > "05ASSEMBLY/ORFs_nucl/02card/${SAMPLE}_ARGs.fa"
done

# Combine all ARCs into a single dataset
mkdir -p 07ARCs_analysis
cat 05ASSEMBLY/ORFs_nucl/01SARG/*_ARGs.fa > 07ARCs_analysis/DWTP_SARG_ARCs.fa
cat 05ASSEMBLY/ORFs_nucl/02card/*_ARGs.fa > 07ARCs_analysis/DWTP_card_ARCs.fa
cat 07ARCs_analysis/DWTP_SARG_ARCs.fa 07ARCs_analysis/DWTP_card_ARCs.fa > 07ARCs_analysis/DWTP_ARCs.fasta

echo "[INFO] ARC extraction completed."