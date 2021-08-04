chrom=6
/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/programs/gcta_1.93.2beta/gcta64 \
--bfile /data/scratch/DGE/DUDGE/MOPOPGEN/plaw/CRC_GWAS/oncoarray/imputation/imputed_plink/CRC_ONCO_chr"$chrom"_imputed \
--chr "$chrom" \
--maf 0.005 \
--cojo-file /data/scratch/DGE/DUDGE/MOPOPGEN/plaw/CRC_GWAS/US_meta/cojo//input/US_UK_ASN_CRC_info4_mega_merged_chr"$chrom".ma \
--cojo-cond /data/scratch/DGE/DUDGE/MOPOPGEN/plaw/CRC_GWAS/US_meta/cojo//new_snps_20200616_"$chrom"_redo.txt \
--out /data/scratch/DGE/DUDGE/MOPOPGEN/plaw/CRC_GWAS/US_meta/cojo//output/EUR/US_UK_ASN_CRC_chr"$chrom"_new_snps_20200616_conditional_redo


