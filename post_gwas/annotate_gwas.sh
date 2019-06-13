#!/bin/bash
#BSUB -J annotate_gwas
#BSUB -P DMPGJMAAQ
#BSUB -W 168:00
#BSUB -o annotate_%J_stdout.txt
#BSUB -e annotate_%J_stderr.txt

#module loads
module load bedtools
module load plink
module load python/3.5.1

outfile="crc_us_meta_info8_0.05"
outpath="/scratch/DGE/MOPOPGEN/plaw/CRC_GWAS/US_meta/meta_out"
#up to the chr mark
meta_data="US_UK_CRC_meta_"
#P-value to filter the meta data on
p_threshold=0.05

cd $outpath

#1) Extract significant results (change p-value to whatever you fancy). If you want to count the number of supporting SNPs in LD then a lower significance value is better.
for chr in {1..22}; do 
	zcat $meta_data$chr\.txt.gz | awk -v pthresh=$p_threshold '{if ($6<pthresh) print $0}' 
done > $outfile\.txt
#get the header from chromosome 1
zcat $meta_data\1.txt.gz | head -1 > head.txt

#2)  Create LD clump input file 
echo "SNP P" > $outfile\_input_plink.txt
cut -d " " -f 2,6 $outfile\.txt >> $outfile\_input_plink.txt

mkdir -p plink_res

#3) Calculate plink clumps within significant hits:
for i in {1..22}; do 
	plink --bfile "/scratch/DGE/MOPOPGEN/plaw/reference_data/1kg_uk10k_eur/1kg_uk10k_eur_chr"$i"_with_cms" --clump $outfile\_input_plink.txt --clump-kb 500 --clump-r2 0.1 --clump-allow-overlap --out plink_res/$outfile\_clump_$i
done

cat plink_res/$outfile\_clump_*.clumped > plink_res/$outfile\_clump_all.clumped

#4) LD WITH GWAS CATALOG (you don’t need to re-run this unless you change the list of known snps)
#for i in `seq 1 22`; do plink --bfile "/scratch/DGE/MOPOPGEN/plaw/reference_data/1kg_uk10k_eur/1kg_uk10k_eur_chr"$i"_with_cms" --r2 --ld-snp-list /data/cancgene-guard2/sdobbins/CRC_GWAS_Jan_2018/LD/GWAS_Catalog_CRC_rsids.txt --out plink_res/LD_$i; done
 
#5) LD WITH KNOWN SNPS (you don’t need to re-run this unless you change the list of known snps)
#ld snps from all published studies - /scratch/DGE/MOPOPGEN/plaw/CRC_GWAS/scotland_meta/output/crc_published_snps.txt
for i in `seq 1 22`; do plink --bfile "/scratch/DGE/MOPOPGEN/plaw/reference_data/1kg_uk10k_eur/1kg_uk10k_eur_chr"$i"_with_cms" --r2 dprime --ld-snp-list /scratch/DGE/MOPOPGEN/plaw/CRC_GWAS/US_meta/known_snps_eur_plink.txt --ld-window 100000 --ld-window-r2 0.01 --ld-window-kb 1000 --out plink_res/LD_known_$i; done
cat plink_res/LD_*.ld | grep -v CHR_A > plink_res/LD_all.ld

#6) Combine data, make an output file with all the annotations
python /scratch/DGE/MOPOPGEN/plaw/scripts/crc/post_gwas/annotate_gwas_results.py $outfile $outpath

#7) Cleanup - remove temporary files
rm head.txt
rm plink_res/$outfile\_clump_all*
rm $outfile\_input_plink.txt
rm $outfile\.txt
