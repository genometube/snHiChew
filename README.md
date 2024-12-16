# snHiChew: High efficient chromatin conformation capture with post-enrichment in single cells

## Publication
Chen Z, Xie Y, Tan C, Ruan F, Luo M, Zhang C, Guo M, Fang Y, Tang C. High efficient chromatin conformation capture without pre-enrichment (HiChew) in single cells. bioRxiv. 2024:2024-06. https://www.biorxiv.org/content/10.1101/2024.06.25.600609v1.abstract

\. Source data:
\. Processed data:


## Overview
snHiChew is a noval single-cell sequencing method that combines efficient sticky-end ligation with post-PCR enrichment using methylation-based selection, addressing key limitations of existing chromatin conformation capture methods. snHiChew demonstrates superior performance with 45-50% valid pair ratios and the ability to generate up to 7.3 million unique valid contacts per cell. The method's scalability and cost-effectiveness make it particularly suitable for large-scale single-cell chromatin conformation studies.

This repository contains all the code that was used to preprocess and analyze snHiChew data from raw fastq. At the start of each file, dependencies like essential libraries or fixed paths to software and data files are declared, which are required for the scripts to function properly.

## Workflow
This workflow utilizes bash scripts to perform snHiChew data analysis, incorporating tools such as HiC-Pro, Cooltools, HiCExplorer, and Higashi.

Demultiplexing raw reads to each barcode is performed by in-house Python scripts. The multi-thread functions are adapted from scHicDemultiplex.py in scHiCExplorer. Must be running on a high I/O speed SSD for optimal performance.

<img src="https://github.com/genometube/snHiChew/blob/main/img/snHiChew.png?raw=true" width="600" height="250">

## Features
Visualization of TAD melting at single cell 50 kb resolution without signal imputation.

<img src="https://github.com/genometube/snHiChew/blob/main/img/melting.png?raw=true" width="600" height="400">

Source data:
NCBI bioproject PRJNA1109567, https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1109567

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
For any questions regarding installation or usage, please contact yemingxie@gmail.com or raise an issue in the GitHub repository.
