#!/bin/bash
#SBATCH --job-name=hic_prep
#SBATCH --cpus-per-task=10

source deactivate
conda activate tools

mkdir -p 01Reads_fq

echo "[INFO] Organizing QC'd Hi-C reads..."
for i in 01READ_QC/*; do 
    if [ -d "$i" ]; then
        b=${i#*/}
        mv "${i}/final_pure_reads_1.fastq" "01Reads_fq/${b}_1.fastq"
        mv "${i}/final_pure_reads_2.fastq" "01Reads_fq/${b}_2.fastq"
    fi
done

echo "[INFO] Counting Hi-C reads..."
for file in 01Reads_fq/*.fastq; do
    echo "$(basename $file .fastq): $(grep -c "^@" $file) reads" >> 01Reads_fq/Hi_C_reads_number.txt
done
echo "[INFO] Preparation completed."