SCRIPTS=/pub/jaclynb1/ERCC1_SingleCell/scripts/slurm_scripts
FPATH1=/pub/jaclynb1/ERCC1_SingleCell/samples22022877
FPATH2=/pub/jaclynb1/ERCC1_SingleCell/samples22032723
OUTPATH=/pub/jaclynb1/ERCC1_SingleCell/Fastqc
FASTQC_DONE=${OUTPATH}/fastqc.done

for sample in ${FPATH1}/*_R*; do
  sbatch ${SCRIPTS}/Step01_slurm_runFastqc.sh $OUTPATH $sample
  echo -e $sample >> $FASTQC_DONE
done

for sample in ${FPATH2}/*_R*; do
  sbatch ${SCRIPTS}/Step01_slurm_runFastqc.sh $OUTPATH $sample
  echo -e $sample >> $FASTQC_DONE
done
