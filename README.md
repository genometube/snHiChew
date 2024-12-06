# snHiChew: High efficient chromatin conformation capture with post-enrichment in single cells

## Overview
snHiChew is a noval single-cell sequencing method that combines efficient sticky-end ligation with post-PCR enrichment using methylation-based selection, addressing key limitations of existing chromatin conformation capture methods. snHiChew demonstrates superior performance with 45-50% valid pair ratios and the ability to generate up to 7.3 million unique valid contacts per cell. The method's scalability and cost-effectiveness make it particularly suitable for large-scale single-cell chromatin conformation studies.

This repository contains all the code that was used to preprocess and analyze snHiChew data from raw fastq to the tabular data. At the start of each file, dependencies like essential libraries or fixed paths to software and data files are declared, which are required for the scripts to function properly.

## Workflow
This workflow performs snHiChew data processing busing bash scripts with HiC-Pro, Cooltools, HiCExplorer, and Higashi encapsulated. 
Demultiplexing is performed by n-house Python scripts. The multi-thread functions are adapted from scHicDemultiplex.py in scHiCExplorer. Must be running on a high I/O speed SSD for optimal performance.

<img src="https://github.com/genometube/snHiChew/blob/main/snHiChew.png?raw=true" width="600" height="260">

## Features
Visualization of TAD melting at single cell 50 kb resolution without signal imputation.

<img src="https://github.com/genometube/snHiChew/blob/main/melting.png?raw=true" width="300" height="400">

## Publication
Processed data:
Source data:

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
