# snHiChew: High efficient chromatin conformation capture with post-enrichment in single cells

## Citations
Chen Z, Xie Y, Tan C, Ruan F, Luo M, Zhang C, Guo M, Fang Y, Tang C. High efficient chromatin conformation capture without pre-enrichment (HiChew) in single cells. bioRxiv. 2024:2024-06. https://www.biorxiv.org/content/10.1101/2024.06.25.600609v1.abstract
All source data can be found in BioProject: PRJNA1109567

● low input bulkHiChew HEK293T 
| Source data | Processed data | 
| --- | --- |
| NCBI biosample SAMN45546211 https://www.ncbi.nlm.nih.gov/biosample/?term=SAMN45546211 | https://drive.google.com/drive/folders/1ep0s-VMBCR0po_w819astXyfCH7SoXje?usp=sharing |

● snHiChew HEK293T 200 cells
| Source data | Processed data | 
| --- | --- |
| NCBI SRA SRR31657844 https://www.ncbi.nlm.nih.gov/sra/?term=SRR31657844 | https://drive.google.com/drive/folders/1AO8DUoMXpcA5-fv7VxKbSw4E3CE306xR?usp=sharing |

● snHiChew GM12878 400 cells
| Source data | Processed data | 
| --- | --- |
| NCBI SRA SRR31657843 https://www.ncbi.nlm.nih.gov/sra/?term=SRR31657843 | https://drive.google.com/drive/folders/1JCeRVXQKgIOiJwU_wucWZZ9QQ4My5HZn?usp=drive_link |

● snHiChew HEK293T wild type 1200 cells
| Source data | Processed data | 
| --- | --- |
| NCBI biosample SAMN41378974 https://www.ncbi.nlm.nih.gov/biosample/?term=SAMN41378974 | https://drive.google.com/drive/folders/16J0pTgnh8IMSAIWTiFffMDWsh09O29fB?usp=sharing |

● snHiChew HEK293T CTCF knockdown 1200 cells
| Source data | Processed data | 
| --- | --- |
| NCBI biosample SAMN42018263	https://www.ncbi.nlm.nih.gov/biosample/?term=SAMN42018263 | https://drive.google.com/drive/folders/1Me_xKag3Ddx8yy8_hhINf4RAJK3FjShb?usp=drive_link |

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
