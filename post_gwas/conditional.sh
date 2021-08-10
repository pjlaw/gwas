chrom=6
#run the conditional analysis
/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/programs/gcta_1.93.2beta/gcta64 \
--bfile /data/scratch/DGE/DUDGE/MOPOPGEN/plaw/CRC_GWAS/oncoarray/imputation/imputed_plink/CRC_ONCO_chr"$chrom"_imputed \
--chr "$chrom" \
--maf 0.005 \
--cojo-file gwas_chr"$chrom".ma \
--cojo-cond snplist_"$chrom".txt \
--out gwas_conditional_chr"$chrom"

conda activate plink1.90b6.18
#extract the conditioned p-values
echo "SNP P" > gwas_conditional_chr"$chrom"_plink.in
awk '{if (NR>1) print $2,$13}'  gwas_conditional_chr"$chrom".cma.cojo >> gwas_conditional_chr"$chrom"_plink.in
#clump the pvalues using plink - keep only GWS snps
plink --bfile /scratch/DGE/MOPOPGEN/plaw/reference_data/1kg_uk10k_eur/1kg_uk10k_eur_chr"$chrom"_with_cms --clump gwas_conditional_chr"$chrom"_plink.in --clump-p1 5e-8 --clump-p2 5e-8  --clump-kb 1000 --clump-r2 0.1 --out gwas_conditional_clump_chr"$chrom"


