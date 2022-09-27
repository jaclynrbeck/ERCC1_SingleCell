#!/bin/bash
#SBATCH --job-name=JB_CR-%J  ## Name of the job.
#SBATCH -p free		     ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=16   ## number of cores the job needs
#SBATCH --error=slurm-%J.err ## error log file
#SBATCH -t 01-00:00:00       ## 1-day run time limit
#SBATCH --mail-type=fail,end
#SBATCH --mail-user=jaclynb1@uci.edu

export PATH=/pub/jaclynb1/cellranger-6.0.1:$PATH
OUTPATH=/pub/jaclynb1/ERCC1_SingleCell/CellRanger
cd $OUTPATH
CR_DONE=${OUTPATH}/cellranger.done

cellranger count --id=${1} --sample=${2} --fastqs=${3} --transcriptome=${4} \
                 --include-introns --expect-cells=10000 --localcores=16 --localmem=32 --no-bam
echo -e $1 >> $CR_DONE


