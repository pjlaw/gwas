#!/bin/bash
# BSUB -P DMPGJMAAQ
# BSUB -J phasing[1-22]
# BSUB -n 3
# BSUB -R "span[hosts=1]"
# BSUB -W 160:00
# BSUB -o /scratch/DGE/MOPOPGEN/plaw/PRS/CRC/onco_with_rel/impute/phasing_%J_chr%I.stdout
# BSUB -e /scratch/DGE/MOPOPGEN/plaw/PRS/CRC/onco_with_rel/impute/phasing_%J_chr%I.stderr


chr=${LSB_JOBINDEX} 

#This script uses shapeit3, which is much faster
#If you'd prefer to use shapeit 2, uncomment the module load, and change command to shapeit


#module load shapeit

/scratch/DGE/MOPOPGEN/plaw/programs/shapeit3/shapeit3.r884.1 -B /scratch/DGE/MOPOPGEN/plaw/PRS/CRC/onco_with_rel/impute/prephased/crc_onco_relatives_chr$chr -M /scratch/Dalek/cancgene/mhenrion/1kG_Dec2013_integrated_unzipped/genetic_map_chr$chr\_combined_b37.txt -O /scratch/DGE/MOPOPGEN/plaw/PRS/CRC/onco_with_rel/impute/phased/crc_onco_relatives_chr$chr\.phased -L /scratch/DGE/MOPOPGEN/plaw/PRS/CRC/onco_with_rel/impute/scripts/tmp/phasing_chr$chr --thread 3
