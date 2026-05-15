# AMR-Phage-Dynamics-AD
**Bioinformatics pipeline for deciphering the dual role of Antiviral Defense Systems (ADS) in shaping environmental resistomes under phage pressure.**

This repository contains the complete custom scripts, data processing pipelines, and methodological configurations supporting our *iMeta* manuscript. It demonstrates how we integrated short/long-read metagenomics, meta-transcriptomics, and high-throughput chromosome conformation capture (Hi-C) to track ARG dynamics in Anaerobic Digestion (AD) systems.

## 🚀 Pipeline Modules

The analytical workflow is modularized into 8 consecutive steps. Please navigate to each directory for detailed scripts, parameter configurations, and methodological rationale.

1. [**Step 1: Quality Control, Assembly, and Binning**](./01_QC_Assembly_Binning)
2. [**Step 2: ARG Quantification, 16S Normalization, and Taxonomic Profiling**](./02_ARGs_Quantification)
3. [**Step 3: Contig-Level Annotation, Mobility Assessment, and Taxonomic Classification**](./03_Gene_Annotation_Mobility)
4. [**Step 4: Viral Sequence Identification, Clustering, Taxonomy, and Novelty Assessment**](./04_Viral_Discovery_and_Taxonomy)
5. [**Step 5: RNA Virus Identification, Clustering, Taxonomy, and Host Prediction**](./05_RNA_Viral_Discovery)
6. [**Step 6: Hi-C Analysis for On-going Phage-Host Interactions**](./06_Hi_C_Interactions)
7. [**Step 7: *In-Memory* Phage-Host Linkage Analysis (DNA Viruses)**](./07_Phage_Host_Linkage_InMemory)
8. [**Step 8: MAG De-replication, Taxonomy, Functional Annotation, and ADS Profiling**](./08_MAGs_Annotation_and_ADS)

## 💡 Note on Reproducibility
To maximize the recovery of viral "dark matter" and robustly define phage-host interactions, this pipeline employs specific union-based ensemble strategies, strict MIUViG clustering thresholds, and tri-evidence host predictions. Detailed parameter choices and their ecological rationales are provided within the `README.md` of each respective module.
