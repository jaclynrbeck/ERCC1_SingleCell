#!/bin/bash

SCRIPTS=/pub/jaclynb1/ERCC1_SingleCell/scripts

sbatch ${SCRIPTS}/slurm_scripts/Step03_slurm_runCellRangerAggr.sh ERCC1Aggr ${SCRIPTS}/cellranger_aggr.csv
