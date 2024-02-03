srun -t 24:00:00 --nodes=1 --cpus-per-task=24 --mem=100G --pty /bin/bash
mkdir /scratch/j.selwyn/asv_classification
cd /scratch/j.selwyn/asv_classification

#1 - make database
module load singularity
BLCA=/work/vollmer/software/blca.sif

singularity exec --bind /work,/scratch,/tmp ${BLCA} \
  1.subset_db_acc.py

#2 - split fasta into NCPU files
#singularity exec --bind /work,/scratch,/tmp ${BLCA} \
#  pyfasta split -n ${SLURM_CPUS_PER_TASK} all_asvs.fasta
#breaks - try without and hope

#3 - Classify ASVs
singularity exec --bind /work,/scratch,/tmp ${BLCA} \
  2.blca_main.py \
  -i all_asvs.fasta \
  -p ${SLURM_CPUS_PER_TASK}
