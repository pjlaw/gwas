#!/bin/bash
# BSUB -P prepay-houlston
# BSUB -J combine[1-22]
# BSUB -n 1
# BSUB -W 100:00
# BSUB -o /scratch/cancgene/plaw/scripts/tmp/combined_%J_chr%I.stdout
# BSUB -e /scratch/cancgene/plaw/scripts/tmp/combined_%J_chr%I.stderr

chr=${LSB_JOBINDEX}

outpath=/scratch/cancgene/plaw/CLL_newdata
mkdir -p $outpath\/combined/

outfile=$outpath\/combined/CLL_newdata_chr$chr\_imputed.gen.gz

if [ -f $outfile ]; then rm $outfile; fi    

#sort by the first chunk number - edit if necessary (split by _ character)!
for file in `ls $outpath/imputed/chr$chr/Main/*.txt | sort -k6,6g -t '_'`; do
	cat $file | gzip >> $outfile
	echo $file
done

echo "";
	
