# Step 1: Quality Control, Assembly, and Binning

This module contains the bioinformatic pipelines for processing raw sequencing data, including Illumina short-read metagenomics/meta-transcriptomics and Oxford Nanopore Technologies (ONT) long-reads. 

## 📁 Directory Structure
Ensure your working directory is structured as follows before running the scripts:
```text
.
├── 01RAW_READS/         # Raw Illumina fastq.gz and ONT fastq files
├── 01READ_QC/           # Output directory for quality-controlled reads
├── 01Reads_fq/          # Centralized directory for clean reads
├── 05ASSEMBLY/          # Metagenomic assembly outputs
└── 06MAGs/              # Binning and refined MAGs outputs

# Step 2: ARG Quantification, 16S Normalization, and Taxonomic Profiling

This module details the read-level bioinformatics analysis. It involves aligning clean reads against the CARD database to profile Antimicrobial Resistance Genes (ARGs), extracting 16S rRNA reads via Kraken2 for subsequent normalization (eliminating depth and biomass biases), and determining high-resolution community taxonomy using MetaPhlAn 4.

## 📁 Directory Structure
```text
.
├── 01Reads_fq/                 # Clean reads from Step 1
├── 04reads_based_analysis/     
│   ├── 02card/                 # CARD alignment outputs and count tables
│   └── 07kraken2/              # Kraken2 16S outputs
└── 10metaphlan4/               # MetaPhlAn4 taxonomic profiles and diversity metrics

# Step 3: Contig-Level Annotation, Mobility Assessment, and Taxonomic Classification

This module performs deep profiling at the contig level. To explicitly link Antimicrobial Resistance Genes (ARGs) with their genomic context (plasmids, ICEs, and host taxonomy), we identified ARG-carrying contigs (ARCs), removed redundancies, and subjected them to a suite of downstream analyses. This step is critical for evaluating the potential for horizontal gene transfer (HGT).

## 📁 Directory Structure
```text
.
├── 05ASSEMBLY/                 # Assembly results from Step 1
│   ├── ORFs_genes/             # Prodigal gene predictions
│   ├── ORFs_protein/           # Prodigal protein sequences (.pro.fa)
│   └── ORFs_nucl/              # Prodigal nucleotide sequences (.fa)
└── 07ARCs_analysis/            # Centralized output for ARC analysis
    ├── 01SARG/                 # SARG annotation
    ├── 02card/                 # CARD annotation
    ├── 04ICEs/                 # ICEs (Integrative & Conjugative Elements)
    ├── 05MGEs/                 # Other MGEs annotation
    ├── 07CAT/                  # Contig Annotation Tool (Taxonomy)
    ├── 08plasflow/             # Plasmid prediction
    └── 11metawrap-classify/    # Taxator-tk taxonomy classification

# Step 4: Viral Sequence Identification, Clustering, Taxonomy, and Novelty Assessment

This module outlines the comprehensive pipeline for identifying viral sequences from both Illumina short-reads and Oxford Nanopore (ONT) long-reads, clustering them into viral Operational Taxonomic Units (vOTUs), and characterizing their taxonomy, lifestyle, and novelty.

> 💡 **Methodological Note :**
> **Identification Strategy (Union vs. Intersection):** To maximize the recovery of viral dark matter in complex anaerobic digestion systems, we adopted a **Union strategy** for initial viral identification using three independent algorithms: VirSorter2 (signature-based), VirFinder (k-mer-based), and geNomad (neural network/marker-based). Requiring an intersection of all three would lead to massive false negatives, especially for novel phages lacking known reference signatures. 
> **Quality Assurance:** To counteract potential false positives introduced by the union approach, we implemented stringent downstream quality controls: 
> 1. Only contigs **≥ 5,000 bp** were considered.
> 2. Viral candidates were clustered at **95% ANI and 85% alignment fraction** (MIUViG guidelines).
> 3. Final validation was performed using **CheckV**, strictly retaining only contigs with identified viral genes or those lacking host genes with lengths >10kb.

## 📁 Directory Structure
```text
.
├── 06VCs/00Contigs/            # Filtered contigs (>5,000 bp)
├── 03_Viral_Discovery/         # Output from VS2, VirFinder, geNomad
├── 01Clustergenome/            # Clustering outputs (vOTUs)
├── 03checkv/                   # Quality control outputs
├── 06ARGs/ & 07VFs/            # ARG and Virulence Factor profiling on vOTUs
└── 04vcontact2/, 08majority/   # Taxonomic classification networks and rules

# Step 5: RNA Virus Identification, Clustering, Taxonomy, and Host Prediction

This module details the pipeline designed specifically for recovering and analyzing RNA viruses from meta-transcriptomic data. Given the distinct biological properties of RNA viruses (e.g., smaller genome sizes and high mutation rates), specific algorithms and modified thresholds were employed compared to the DNA virus pipeline.

> 💡 **Methodological Note (Response to Reviewer Guidelines):**
> * **Length Threshold:** Contigs ≥ 1,000 bp were retained (instead of 5,000 bp) due to the intrinsically smaller genome sizes of RNA viruses.
> * **Identification Strategy:** We coupled a highly conserved signature-based approach (HMM profiles of RNA-directed RNA polymerase, RdRp) with a machine-learning framework (geNomad Riboviria models) to maximize RNA virus recovery while maintaining high precision.
> * **Clustering Parameters:** Following MIUViG guidelines specifically for RNA viruses, vOTUs were clustered using slightly relaxed parameters: **≥90% ANI and ≥80% alignment fraction**, accommodating their higher evolutionary mutation rates.

## 📁 Directory Structure
```text
.
├── 08RNA_VCs/
│   ├── 00Contigs/           # RNA contigs (>1,000 bp) and ORFs
│   ├── 01genomad/           # geNomad prediction outputs
│   ├── rdrp/                # RdRp-search module outputs
│   ├── 01Clustergenome/     # RNA vOTUs clustering outputs
│   └── 10host-link/         # Virus-Host interaction networks (CRISPR, EVEs, RNAVirHost)

# Step 6: Hi-C Analysis for On-going Phage-Host Interactions

This module details the processing of high-throughput chromosome conformation capture (Hi-C) data to establish *in vivo*, active phage-host interactions. By chemically cross-linking DNA physically interacting within intact cells prior to extraction, Hi-C allows us to capture the exact moment a phage genome is inside a bacterial host.

> 💡 **Methodological Note (Response to Reviewer Guidelines):**
> * **Mapping Parameters:** We utilized `bwa mem -5SP`. The `-5SP` flag is critical for Hi-C data as it forces strict paired-end mapping and skips rescuing secondary/supplementary alignments, which could otherwise create artificial chimeric cross-links.
> * **Strict Filtering:** Mapped reads were filtered with `samtools view -F 0x904` to eliminate unmapped reads, secondary alignments, and supplementary alignments.
> * **Thresholding:** To completely eliminate spurious noise and random ligations inherent in bulk Hi-C library preparation, we strictly required a minimum of **>5 distinct Hi-C paired-end links** between a viral contig and a Metagenome-Assembled Genome (MAG/bin) to confirm an *on-going* infection event.

## 📁 Directory Structure
```text
.
├── 01READ_QC/           # Raw QC'd Hi-C reads
├── 01Reads_fq/          # Organized fastq files for mapping
├── 02Hi-C/              # BWA mapping, SAM/BAM conversions, and interaction networks
└── hic_assembly.fasta   # Combined reference database (All MAGs + All vOTUs)

# Step 7: *In-Memory* Phage-Host Linkage Analysis (DNA Viruses)

While Hi-C (Step 6) captures *on-going* physical infections, this module reconstructs the *in-memory* (historical) infection networks between DNA phages and their bacterial hosts. Bacterial genomes record past phage encounters through adaptive immune systems and genetic exchanges. 

> 💡 **Methodological Note (Response to Reviewer Guidelines):**
> We employed a robust **tri-evidence approach** to computationally predict these historical linkages:
> 1. **CRISPR-Cas Spacers:** Extremely high confidence. We extracted spacers from host MAGs and mapped them against viral contigs using strict short-read alignment parameters (`blastn-short`, 97% identity).
> 2. **tRNA Matches:** Phages often hijack or integrate near host tRNAs. We identified exact matches (100% identity) between viral tRNAs and host genomes.
> 3. **Sequence Homology:** Broad homologous regions (≥80% identity) indicating prophage integration, horizontal gene transfer, or recombination events.

## 📁 Directory Structure
```text
.
├── 06MAGs/04dRep/           # Dereplicated high-quality MAGs
├── 01Clustergenome/         # Viral Operational Taxonomic Units (vOTUs)
└── 10host-link/             # Output directories for linkages
    ├── 01crispr/            # CRISPR spacer extraction and mapping
    ├── 02tRNA/              # tRNA extraction and mapping
    └── 03homology/          # Genomic homology mapping

# Step 8: MAG De-replication, Taxonomy, Functional Annotation, and ADS Profiling

This final module operates at the Metagenome-Assembled Genome (MAG) level. It establishes a high-quality, non-redundant catalog of microbial genomes from the anaerobic digestion (AD) system, determines their precise taxonomy, profiles their Antimicrobial Resistance Genes (ARGs) and Mobile Genetic Elements (MGEs), and comprehensively maps their Antiviral Defense Systems (ADS).

> 💡 **Methodological Note (Response to Reviewer Guidelines):**
> * **Dereplication Strategy:** We used `dRep` to dereplicate bins across all samples at a 99% primary and 95% secondary ANI threshold. This ensures we are analyzing unique, species-level genomic representatives rather than redundant strain variations, avoiding statistical bias in host abundance.
> * **ADS Profiling (Dual-Tool Approach):** To rigorously evaluate the hypothesis that ARBs accumulate defense systems to persist under intense phage pressure, we utilized both **PADLOC** and **DefenseFinder**. Employing this dual-tool approach ensures maximum sensitivity and comprehensive classification of the diverse immune repertoires (e.g., CRISPR-Cas, R-M, ABI) encoded by these environmental genomes.

## 📁 Directory Structure
```text
.
├── 06MAGs/
│   ├── 02refine_bins/        # Refined bins from MetaWRAP (Step 1)
│   ├── 04dRep/               # Non-redundant MAGs catalog
│   ├── 06gtdbtk_classify/    # GTDB-Tk taxonomy output
│   ├── ORFs_nucl/            # Functional annotation (ARGs, ICEs, VFs, etc.)
│   └── Defense_Systems/      # PADLOC and DefenseFinder outputs
