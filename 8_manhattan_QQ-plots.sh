#!/bin/bash
#BSUB -n 1
#BSUB -W 120:00
#BSUB -P prepay-houlston
#BSUB -o /scratch/cancgene/mhenrion/RCC_for_phasing_imputing/NCI_GWAS/filtered_results/combined/NCIGWAS_RCC_results_check_stdout.txt
#BSUB -e /scratch/cancgene/mhenrion/RCC_for_phasing_imputing/NCI_GWAS/filtered_results/combined/NCIGWAS_RCC_results_check_stderr.txt

R --vanilla < /scratch/cancgene/mhenrion/RCC_for_phasing_imputing/NCI_GWAS/manhattan_QQ-plots.R --args /scratch/cancgene/mhenrion/RCC_for_phasing_imputing/NCI_GWAS/filtered_results/combined/NCIGWAS_RCC_results.txt /scratch/cancgene/mhenrion/RCC_for_phasing_imputing/NCI_GWAS/filtered_results/combined/NCIGWAS_RCC_results_check


