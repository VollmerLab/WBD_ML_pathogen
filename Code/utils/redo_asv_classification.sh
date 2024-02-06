#Login to Discovery HPC
srun -t 24:00:00 --nodes=1 --cpus-per-task=24 --mem=100G --pty /bin/bash
outdir=/scratch/j.selwyn/asv_classification #Change to where you want to use as a working directory
fasta_file=all_asvs.fasta #FASTA File with all 16s squences

mkdir -p ${outdir}
cd ${outdir}
##put the ${fasta_file} into the ${outdir}

#1 - make database
module load singularity
BLCA=/work/vollmer/software/blca.sif

singularity exec --bind /work,/scratch,/tmp ${BLCA} \
  1.subset_db_acc.py

#2 - Classify ASVs
singularity exec --bind /work,/scratch,/tmp ${BLCA} \
  2.blca_main.py \
  -i ${fasta_file} \
  -p ${SLURM_CPUS_PER_TASK}
