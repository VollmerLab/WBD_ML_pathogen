Get other WBD 16s sequencing
```
mkdir -p /scratch/j.selwyn/old_wbd
cd /scratch/j.selwyn/old_wbd
```

## Download fastqs
(Gignoux-Wolfsohn et al. 2017) - PRJNA387312 - Panama July 2014
`bash download_gignoux-wolfohn2017.sh`

(Gignoux-Wolfsohn and Vollmer 2015) - SRR2078055-SRR2078108 - Panama 2009/2010
```
R
write_lines(str_c('SRR', 2078055:2078108),
            '../../intermediate_files/gw_2015.txt')
```

(Rosales et al. 2019) - Acerv/Apalm - Florida 2017
- https://figshare.com/articles/dataset/Untitled_Item/8226209?file=15330449

(Gignoux-Wolfsohn et al. 2020) - PRJNA511881 - acerv - florida 2013/2014 october


(Certner and Vollmer 2018) -
``

(Pantos & Bythell 2006) - AY323132–AY323197 - Apalm - Barbados
```
R
write_lines(str_c('AY', 323132:323197),
            '../../intermediate_files/pantos_bythell_2006.txt')
```
Copy to https://www.ncbi.nlm.nih.gov/sites/batchentrez
Download and rename to `pantos_2006.fasta`

## Merge fastqs
gw_2017
```

```

## Copy 16s database
```
cp -r /scratch/j.selwyn/asv_classification/db /scratch/j.selwyn/old_wbd
```
## Classify Again
```
srun -t 24:00:00 --nodes=1 --cpus-per-task=24 --mem=100G --pty /bin/bash
cd /scratch/j.selwyn/old_wbd
module load singularity
BLCA=/work/vollmer/software/blca.sif

singularity exec --bind /work,/scratch,/tmp ${BLCA} \
  2.blca_main.py \
  -i pantos_2006.fasta \
  -p ${SLURM_CPUS_PER_TASK}

singularity exec --bind /work,/scratch,/tmp ${BLCA} \
  2.blca_main.py \
  -i rosales_2019.fasta \
  -p ${SLURM_CPUS_PER_TASK}

singularity exec --bind /work,/scratch,/tmp ${BLCA} \
  2.blca_main.py \
  -i klinges_2022.fasta \
  -p ${SLURM_CPUS_PER_TASK}

singularity exec --bind /work,/scratch,/tmp ${BLCA} \
  2.blca_main.py \
  -i trytten_unpub.fasta \
  -p ${SLURM_CPUS_PER_TASK}
```
