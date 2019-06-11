#!/bin/bash
#BSUB -J filter_snptest[1-22]
#BSUB -P prepay-houlston
#BSUB -n 1
#BSUB -W 120:00
#BSUB -o /scratch/cancgene/plaw/scripts/tmp/filter_%J.chr%I_stdout.txt
#BSUB -e /scratch/cancgene/plaw/scripts/tmp/filter_%J.chr%I_stderr.txt

# set-up run parameters
chr=${LSB_JOBINDEX}

mafFilter=0.005
infoFilter=0.4
hweFilter=1e-5

path=/scratch/cancgene/plaw/CLL_newdata/
#mkdir -p $path\/filtered/
inFile=$path\/snptest/CLL_newdata_snptest_chr$chr\.txt.gz 
outFile=$path\/filtered/CLL_newdata_snptest_filtered_chr$chr\.txt 
lineSkip=13
zcat $inFile | head -$lineSkip | tail -1 > $outFile
#filter cases and controls maf separately, also filter both on general info and model info score
gzip -dc $inFile | awk -v hweFilter=$hweFilter -v mafFilter=$mafFilter -v infoFilter=$infoFilter -v lineSkip=$lineSkip '{if(FNR>lineSkip && $9>infoFilter && $46>infoFilter && $30>mafFilter && $31>mafFilter && $35>hweFilter && $45!="NA") print }' >> $outFile 
gzip $outFile
	
