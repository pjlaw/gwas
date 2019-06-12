#!/bin/bash
#BSUB -J pca
#BSUB -P DMPGJMAAQ
#BSUB -q normal
#BSUB -R "span[hosts=1]"
#BSUB -n 2
#BSUB -W 160:00
#BSUB -o /scratch/DGE/MOPOPGEN/plaw/scripts/tmp/pca_%J_stdout.txt
#BSUB -e /scratch/DGE/MOPOPGEN/plaw/scripts/tmp/pca_%J_stderr.txt

module load eigensoft

ipath="/scratch/DGE/MOPOPGEN/plaw/CRC_GWAS/oncoarray/imputation/unphased/redoQC"
study="crc_onco_caco_clean_hapmap_merged"

prefix=$ipath\/$study
#~ cp $prefix\.bim $prefix\.pedsnp
#~ cp $prefix\.fam $prefix\.pedind

perl /scratch/DGE/MOPOPGEN/plaw/programs/smartpca_v6.pl -i $prefix\.bed -a $prefix\.pedsnp -b $prefix\.pedind -o $prefix\.pca -p $prefix\.pca.plot -e $prefix\.pca.eval -l $prefix\.pca.log -k 3 -t 3 -w /scratch/DGE/MOPOPGEN/plaw/Dalek/scripts/qc/pca-populations.txt
