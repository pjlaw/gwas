#!/bin/bash
#SBATCH -J meta
#SBATCH -t 120:00:00
#SBATCH --array=1-22
#SBATCH -o meta.%A_%a.out

chr=${SLURM_ARRAY_TASK_ID} 

path="/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/data_in/"
/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/programs/meta_v1.7_x86_64_dynamic/meta --method 1 --threshold 0.8 --cohort $path/CLL1_58BC/filtered/CLL1_filtered_allcols_all_maf_chr$i\.txt $path/CLL2_NBS/filtered/CLL2_filtered_allcols_chr$i\.txt --output /scratch/cancgene/plaw/CLLmeta/meta_random/CLL_meta_random_chr$i\.txt
