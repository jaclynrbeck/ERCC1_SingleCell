#!/bin/bash
#SBATCH --job-name=JB_Fastqc      ## Name of the job.
#SBATCH -p free		          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=1    ## number of cores the job needs
#SBATCH --error=slurm-%J.err ## error log file
#SBATCH --mail-type=fail,end
#SBATCH --mail-user=jaclynb1@uci.edu

module load fastqc/0.11.9

fastqc --noextract --threads 1 --outdir ${1} ${2}
