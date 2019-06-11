#!/bin/bash
#BSUB -J CLL_meta[1-22]
#BSUB -P zero-houlston
#BSUB -n 1
#BSUB -W 120:00
#BSUB -o /scratch/cancgene/plaw/scripts/tmp/CLL_%J.chr%I_meta_stdout.txt
#BSUB -e /scratch/cancgene/plaw/scripts/tmp/CLL_%J.chr%I_meta_stderr.txt

# set-up run parameters
i=${LSB_JOBINDEX}
path="/data/sutst/DGE/MOPOPGEN/DATA/studies/CLL/impute_UK10k_1kG"

/scratch/cancgene/mhenrion/meta_snptest/meta_static/meta-stat-1.3.2 --snptest --method 1 --threshold 0.4 --cohort $path/CLL1_58BC/filtered/CLL1_filtered_allcols_all_maf_chr$i\.txt $path/CLL2_NBS/filtered/CLL2_filtered_allcols_chr$i\.txt --output /scratch/cancgene/plaw/CLLmeta/meta_random/CLL_meta_random_chr$i\.txt
