#!/bin/bash

source deactivate
conda activate tools

# Process each ONT fastq file
for F in 01RAW_READS/*.fastq; do
    BASE=${F##*/}    
    SAMPLE=${BASE%.*}
    
    mkdir -p 01READ_QC/$SAMPLE   
    # Run fastplong for long-read quality control
    fastplong -i $F -o 01READ_QC/$SAMPLE/"$SAMPLE".fastq
    
    # Organize output reports and clean sequences
    mv fastplong.html 01READ_QC/$SAMPLE/
    mv fastplong.json 01READ_QC/$SAMPLE/
    mv 01READ_QC/$SAMPLE/"$SAMPLE".fastq 01Reads_fq/
done