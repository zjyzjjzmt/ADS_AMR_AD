#!/bin/bash
#SBATCH --job-name=viral_discovery_illumina
#SBATCH --cpus-per-task=50

source deactivate
conda activate tools

# 1. Filter contigs > 5000bp
mkdir -p 06VCs/00Contigs
for F in 01Reads_fq/*_1.fastq.gz; do
    BASE=${F##*/}
    SAMPLE=${BASE%_*}
    # Rename contig headers with sample names and filter length
    sed "s/k141/$SAMPLE/g" 05ASSEMBLY/$SAMPLE/final_assembly.fasta | seqtk seq -L 5000 > 06VCs/00Contigs/$SAMPLE.fa
done
cat 06VCs/00Contigs/*.fa > 06VCs/00Contigs/DWTP_contigs5000bp.fasta

CONTIGS_ALL="06VCs/00Contigs/DWTP_contigs5000bp.fasta"
OUT_DIR="./03_Viral_Discovery"
mkdir -p ${OUT_DIR}

# 2. Run VirSorter2 (Parameter: --include-groups dsDNAphage,ssDNA)
virsorter run -i ${CONTIGS_ALL} -w ${OUT_DIR}/vs2_out --include-groups dsDNAphage,ssDNA -j 50 all

# 3. Run geNomad
genomad end-to-end --cleanup ${CONTIGS_ALL} ${OUT_DIR}/genomad_out /path/to/genomad_db

# 4. Run VirFinder
Rscript scripts/run_virfinder.R ${CONTIGS_ALL} ${OUT_DIR}/virfinder_out.tsv

# 5. Extract Union of Candidates and Remove Redundancy
# filter_viral_contigs.py applies the UNION logic (VS2 OR VF OR geNomad)
python scripts/filter_viral_contigs.py \
    --vs2 ${OUT_DIR}/vs2_out/final-viral-score.tsv \
    --vf ${OUT_DIR}/virfinder_out.tsv \
    --genomad ${OUT_DIR}/genomad_out/*_summary.tsv \
    --fasta ${CONTIGS_ALL} \
    --out ${OUT_DIR}/merged_viral_candidates.fasta

cd-hit-est -i ${OUT_DIR}/merged_viral_candidates.fasta -o ${OUT_DIR}/illumina_viral_nr.fasta -c 1.0 -n 10 -M 16000 -T 50