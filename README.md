# ERCC1_SingleCell

This repository contains the code used to process single cell sequencing data 
for my ERCC1 KO thesis project. It re-uses a lot of code and processes from
https://github.com/jaclynrbeck/TCellsAD2022 for consistency. 

Description of folders:

- `cellranger_scripts` - Shell scripts used to run CellRanger on UCI's High 
Performance Computing (HPC) cluster, using Slurm for job management. 
- `data` - Contains count matrices output by CellRanger and several 
publicly-available DEG data sets related to aging and AD from other papers[^bignote]. 
This folder is also where data from the pipeline is output by default. 
- `functions` - Helper functions for various steps in the pipeline. These files 
are copied from https://github.com/jaclynrbeck/TCellsAD2022, including functions
that may not be used in this pipeline. 

---

## (Optional) Run CellRanger

This step has already been performed and the output is in this repository, but
in case CellRanger ever needs to be re-run on the raw data again, follow these
steps:

1. Request the fastq files. 
2. Download the human reference data for gene expression ([refdata-gex-GRCh38-2020-A.tar.gz](https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz)) from 10X genomics. (Optional) Consider trying the [combined mouse + human reference](https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-and-mm10-2020-A.tar.gz). I've tried this and it did label some cells as "mouse" and some
   as doublets, but I haven't looked at that data in Seurat to see if it makes
   sense.
3. Every csv file in `cellranger_scripts` needs to be edited so that file paths
point to your file structure. 
4. Shell scripts will need to be tweaked if not using Slurm or if using a 
different way of configuring Slurm than UCI's HPC. 
5. Run the shell scripts in order from Step01 to Step03.

---

## Run Pipeline

Most of the filenames in `Filenames.R` should be set to good default values, and
file paths are relative to the main directory. Edit filenames as necessary.

Scripts are intended to be run in order from Step01 through Step05. See
comments at the header of individual files for more information about each 
script. Steps 04 and 05 were not used in my thesis but may have some interesting
information in them, so I included them.

Output from the scripts will be put in `data` under various sub-folders.

If run as-is, this pipeline will output the data and figures used for
my thesis. 

[^bignote]: Papers used for reference data sets:  
Galatro 2017: https://doi.org/10.1038/nn.4597  
Holtman 2015: https://doi.org/10.1186/s40478-015-0203-5  
Keren-Shaul 2017: http://dx.doi.org/10.1016/j.cell.2017.05.018  
Lopes 2022: https://doi.org/10.1038/s41588-021-00976-y  
Olah 2018: https://doi.org/10.1038/s41467-018-02926-5  
Srinivasan 2020: https://doi.org/10.1016/j.celrep.2020.107843  
CellAge Database: https://genomics.senescence.info/cells/  
GenAge Database: https://genomics.senescence.info/download.html  

