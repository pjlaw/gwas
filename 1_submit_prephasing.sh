#!/bin/bash

module load shapeit

for chr in {1..22}; do
    echo '#!/bin/bash
    #BSUB -J "shapeit_p2"
	#BSUB -n 1
	#BSUB -P zero-houlston
	#BSUB -W 100:00
	#BSUB -o tmp/output.%J
	#BSUB -e tmp/errors.%J' > tmp/tmp_shapeit_cmd_$chr\.sh
    
    echo "" >> tmp/tmp_shapeit_cmd_$chr\.sh
    echo "shapeit -B /scratch/cancgene/plaw/plink_gwas/unphased/CLL2_chr$chr -M /scratch/cancgene/mhenrion/1kG_Dec2013_integrated_unzipped/genetic_map_chr$chr\_combined_b37.txt -O /scratch/cancgene/plaw/plink_gwas/phased/CLL2_chr$chr\.phased --thread 12" >> tmp/tmp_shapeit_cmd_$chr\.sh
       
    bsub < tmp/tmp_shapeit_cmd_$chr\.sh
    
done;
