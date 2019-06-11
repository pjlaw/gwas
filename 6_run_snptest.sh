#!/bin/bash
# BSUB -P prepay-houlston
# BSUB -J snptest[1-22]
# BSUB -n 1
# BSUB -W 100:00
# BSUB -o /scratch/cancgene/plaw/scripts/tmp/snptest_%J_%I.stdout
# BSUB -e /scratch/cancgene/plaw/scripts/tmp/snptest_%J_%I.stderr

chr=${LSB_JOBINDEX}

path=/scratch/cancgene/plaw/CLL_newdata

impFile=$path\/combined/CLL1_chr$chr\.imputed.gen.gz
sampFile=$path\/CLL1_sample_file.sample.gz
outFile=$path\/snptest/CLL1_snptest_chr$chr\.txt.gz

#run snptest
/scratch/cancgene/plaw/programs/snptest/snptest_v2.5.1_linux_x86_64_static/snptest_v2.5.1 -assume_chromosome $chr -method score -frequentist 1 -pheno pheno_caco -hwe -data $impFile $sampFile -o $outFile 

