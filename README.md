# snHiChew: High efficient chromatin conformation capture with post-enrichment in single cells

## Citations
Chen Z, Xie Y, Tan C, Ruan F, Luo M, Zhang C, Guo M, Fang Y, Tang C. High efficient chromatin conformation capture without pre-enrichment (HiChew) in single cells. bioRxiv. 2024:2024-06. https://www.biorxiv.org/content/10.1101/2024.06.25.600609v1.abstract

All source data can be found in NCBI BioProject: PRJNA1109567

| Key samples | Source data | Processed data | 
| :--- | :--- | :--- |
| bulkHiChew_HEK293T | NCBI biosample SAMN45546211 https://www.ncbi.nlm.nih.gov/biosample/?term=SAMN45546211 | https://drive.google.com/drive/folders/1ep0s-VMBCR0po_w819astXyfCH7SoXje?usp=sharing |
| snHiChew_HEK293T 200_cells | NCBI SRA SRR31657844 https://www.ncbi.nlm.nih.gov/sra/?term=SRR31657844 | https://drive.google.com/drive/folders/1AO8DUoMXpcA5-fv7VxKbSw4E3CE306xR?usp=sharing |
| snHiChew_GM12878 400_cells | NCBI SRA SRR31657843 https://www.ncbi.nlm.nih.gov/sra/?term=SRR31657843 | https://drive.google.com/drive/folders/1JCeRVXQKgIOiJwU_wucWZZ9QQ4My5HZn?usp=drive_link |
| snHiChew_HEK293T_wt 1200_cells | NCBI biosample SAMN41378974 https://www.ncbi.nlm.nih.gov/biosample/?term=SAMN41378974 | https://drive.google.com/drive/folders/16J0pTgnh8IMSAIWTiFffMDWsh09O29fB?usp=sharing |
| snHiChew_HEK293T_CTCFkd 1200_cells | NCBI biosample SAMN42018263	https://www.ncbi.nlm.nih.gov/biosample/?term=SAMN42018263 | https://drive.google.com/drive/folders/1D7tMfyKbqH5SNil0JFXrNGmcHfg1tsB-?usp=sharing |
| ChIP-seq HEK293T_WT_CFCFkd | NCBI biosample SAMN45357173 https://www.ncbi.nlm.nih.gov/biosample/45357173 | https://drive.google.com/drive/folders/1TFbv85nOPkm0HMHTZwsxYOzjexfrTIgz?usp=drive_link |

## Overview
snHiChew is a noval single-cell sequencing method that combines efficient sticky-end ligation with post-PCR enrichment using methylation-based selection, addressing key limitations of existing chromatin conformation capture methods. snHiChew demonstrates superior performance with 45-50% valid pair ratios and the ability to generate up to 7.3 million unique valid contacts per cell. The method's scalability and cost-effectiveness make it particularly suitable for large-scale single-cell chromatin conformation studies.

## System requirements and usage
This repository contains all the code that was used to preprocess and analyze snHiChew data from raw fastq. At the start of each file, dependencies like essential libraries or fixed paths to software (tested version number) and data files are declared, which are required for the scripts to function properly. 

Demultiplexing raw reads to each barcode is performed by in-house Python scripts *demultiplex_hichew_top_cell.py* and *demultiplex_hichew_multicore_gz.py*. The multi-thread functions are adapted from [*scHicDemultiplex.py*](https://github.com/joachimwolff/scHiCExplorer/blob/master/schicexplorer/scHicDemultiplex.py) in scHiCExplorer (v7). Must be running on a high I/O speed SSD for optimal performance.

## Workflow
This workflow utilizes bash scripts to perform snHiChew data analysis, incorporating tools such as HiC-Pro, Cooltools, HiCExplorer, and Higashi.

<img src="https://github.com/genometube/snHiChew/blob/main/img/snHiChew.png?raw=true" width="1000" height="360">

## Instructions for each module

### 01.Demultiplex
Copy all files under [01.Demultiplex](https://github.com/genometube/snHiChew/tree/main/01.Demultiplex) module to your directory for demultiplexing. 

Configure the resources for demultiplexing, loaded cell number, file path, and python path in *demultiplex_i5_i7.sh*. And run *demultiplex_i5_i7.sh*.

---
### 02.HiC-Pro
Copy all files under [02.HiC-Pro](https://github.com/genometube/snHiChew/tree/main/02.HiC-Pro) module to your directory to run HiC-Pro for each cell.

Configure the resources and reference files for HiC-Pro analysis in each cell in *config-hicpro.txt*. To prevent memory errors when processing the full snHiChew dataset, set `JOB_MEM = 20` or higher.

Configure the PATH to demultiplexed single cell fastq and PATH to HiC-Pro in *parallel.sh* and *hic-pro-single.sh*. Then run *parallel.sh*.

---
### 03.QC_metrics
Copy all files under [03.QC_metrics](https://github.com/genometube/snHiChew/tree/main/03.QC_metrics) module to your directory to run QC for HiC-Pro results.

Configure the PATH to bedtools, hic pro reference files, and PATH to hic pro results in *hic-pro-metrics.sh*.

Configure the PATH to demultiplexed single cell fastq in *parallel.sh*. Then run *parallel.sh*.

To obtain valid cells by knee point detection (R package [kneedle](https://github.com/etam4260/kneedle)), run *rm_empty_cell_auto.sh* with arguments that specified sample_id or project_id, knee1_sensitivity, knee2_sensitivity, and knee1_knee2_interval. 

Example command: `sh rm_empty_cell_auto.sh snHiChew_HEKwt 1 2 150`

---
### 04.Pseudo_bulk
Copy all files under [04.Pseudo_bulk](https://github.com/genometube/snHiChew/tree/main/04.Pseudo_bulk) module to your directory.
Configure the compartment bin size and PATH to cooltools in *cooltools-gc.sh*. Then run *cooltools-gc.sh* to create reference-track for cooltools call-compartments.

To generate pooled contacts in `.hic` and `.mcool` format, configure the PATH to softwares dependencies and files listed in *stack_ab_ins_loop.sh*. Then run *stack_ab_ins_loop.sh*.

---
### 05.Melting_analysis
Copy all files under [05.Melting_analysis](https://github.com/genometube/snHiChew/tree/main/05.Melting_analysis) module to your directory.
SubTAD melting was observed when unique validpair is above at least 250k for each cell. To analyze single cell chromatin contact at 50k bin resolution, valid snHiChew cells were first downsampled to the same number of unique valid pairs and converted to mcools using Cooltools to call TAD-separation score (scripts in [validpair_downsample_mcool](https://github.com/genometube/snHiChew/tree/main/05.Melting_analysis/validpair_downsample_mcool)). 

The melting and concretion cell clusters were obtained by clustering on single cell TAD-separation score. Hierarchical clustering (ward.D or ward.D2) was performed on the scaled single-cell TAD-separation score for each 50k genomic bin in the region of interest, which can be found in the plotting_scripts folder [FINAL_ins_downsample_multi_median_split_wardD_20250110.R](https://github.com/genometube/snHiChew/blob/main/ploting_scripts/single_cell_compare/validpair_downsample/cpd/FINAL_ins_downsample_multi_median_split_wardD_20250110.R).

The MELTRON pipeline (PMID: 34789882) was adapted to define the melting state of each genomic region of interest. In brief, melting states were calculated for HEK293T B compartment genes longer than 300 kb on the hg19 genome as the genomic region of interest (scripts in [long_transcript_melting_wardD](https://github.com/genometube/snHiChew/tree/main/05.Melting_analysis/long_transcript_melting_wardD)). The genome wide analysis of melting status for these long genes can be found in the plotting scripts [melting_exp_tx.R](https://github.com/genometube/snHiChew/blob/main/ploting_scripts/long_gene_melting/final_hekWT/melting_exp_tx.R).

---
### 06.Higashi
Codes for higashi analysis can be found under [06.Higashi](https://github.com/genometube/snHiChew/tree/main/06.Higashi). Dimensionality reduction were performed by [Higashi](https://github.com/ma-compbio/Higashi) with standard settings.

To annotation the HEK293T cell phase in snHiChew, the HEK293T 2-phase Repli-seq dataset (4DNESSV33VOL, 4DNESH4XLJCW) was obtained from the 4D Nucleome Data Portal. The early/late repli-score ratio for each valid cell was labeled following the repli-seq methods described in the Nature Protocols (PMID: 29599440) and Nagano et al. 2017 (PMID: 28682332).

---
### Plotting scripts
All plotting codes for the paper can be found in the folder [plotting_script](https://github.com/genometube/snHiChew/tree/main/ploting_scripts).

## Demo
Codes and outputs for snHiChew demo scripts can be found in the folder [demo](https://github.com/genometube/snHiChew/tree/main/demo), which is running on snHiChew GM12878 400 cells datasets, NCBI SRA SRR31657843 https://www.ncbi.nlm.nih.gov/sra/?term=SRR31657843.

## Features
Visualization of TAD melting at single cell 50 kb resolution without signal imputation.

<img src="https://github.com/genometube/snHiChew/blob/main/img/melting.png?raw=true" width="600" height="400">

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
For any questions regarding installation or usage, please contact yemingxie@gmail.com or raise an issue in the GitHub repository.
