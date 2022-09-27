#!/bin/bash

SCRIPTS=/pub/jaclynb1/ERCC1_SingleCell/scripts
CSV=/pub/jaclynb1/ERCC1_SingleCell/scripts/config_cellranger_count_allsamples.csv
REF=/pub/jaclynb1/refs/cellranger/refdata-gex-GRCh38-2020-A

while IFS="," read -r name sample fastq
do
  sbatch ${SCRIPTS}/slurm_scripts/Step02_slurm_runCellRanger.sh ${name} ${sample} ${fastq} ${REF}
done < <(tail -n +2 ${CSV})
