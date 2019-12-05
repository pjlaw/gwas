# GWAS pipeline

Houlston Lab scripts for running a GWAS. Scripts written by Marc Henrion and Philip Law, with contributions from Sara Dobbins, Nicky Whiffin, Yusanne Ma, Yufei Wang

## Quality control
Before you start, first perform quality control on the SNPs. 
We generally follow this protocol:
http://www.nature.com/nprot/journal/v5/n9/abs/nprot.2010.116.html 

There is a newer version of plink which is much faster (it's the default on the HPC). However, some of the commands have changed slightly. 
You can read up on the various commands here: https://www.cog-genomics.org/plink2/ 

Updated scripts in the QC folder.

## 0.0 Remove SNPs that are not in the reference panel
davros:/scratch/DGE/MOPOPGEN/plaw/reference_data/philsrefpanelsnpslist.txt.gz

Generated with */scratch/DGE/MOPOPGEN/plaw/Dalek/scripts/gwas/find_ref_snps_in_data.py* [may need to do this on a highmem node]

Then extract the wanted_snps using plink

Do this in R:
```ref_snps=read.table(gzfile(“philsrefpanelsnpslist.txt.gz”), stringsAsFactors=F)[,1]
geno_snps=read.table(“study.bim”, stringsAsFactors=F)[,2]
wanted_snps=intersect(ref_snps, geno_snps)
cat(wanted_snps, sep=”\n”, file=”wanted_snps.txt”)
```

## 0.1 split your dataset by chromosome using plink
`for chr in {1..22}; do plink --allow-no-sex --bfile YOURFILE --chr $chr --make-bed --out YOUROUTPUTFILE_chr$chr; done`

## 1. phase using shapeit
*1_submit_prephasing.sh*

Change the paths to your data


## 2. impute using imputev2
*2_whole_genome_imputation_UK10K_1kG.pl*

1. create the scripts that will be submitted (this is the updated version of the famous sara-yussanne-yufei-nicky script). You should only need to change the top couple lines (input parameters) 
2. submit these files (sitting in whatever folder you specified as the scriptDir in the above script):
`for file in impute*sh; do bsub < $file; done`

## 3. check the imputation
*3a_check_whole_genome_imputation.pl*

This script will go through the output files and check the imputation finished correctly. It creates an output file (impute_check.log, say) that you can grep to check for error during imputing:
* SUPER (these chunks finished as expected)
* OK (no output files for these chunks, but accounted for in the impute logs; essentially this will be chunks falling entirely within centromeres etc)
* CAUTION (these finished successfully, but the imputation accuracy is not as good as you would have liked)
* WARNING (these are the ones you really do not want -- these finished without output file and no reason given why -- usually what occurs if the scheduler kills your job for exceeding memory etc)
I've also added CAUTION2 and WARNING2 where there was a number of warnings reported.

The parsing done by this script isn't perfect, so if you see something weird, check the actual impute log file.

As a convenience function, you can run
*3a_check_whole_genome_imputation.pl*
which will grep for these categories and split them into separate files (SUPER and OK into the same file)

Re-submit anything that needs to be re-run; if your job exceeded the memory limit, be sure to either lower the number of jobs per node or increase the number of cores (hence memory) requested

## 4. combine per chromosome
*4_combine_imputed_files.sh*

Combine and compress the impute data by chromosome. This script sorts the files by bp order, so make sure you are sorting in the correct field (the -k parameter)

## 5. create a sample file for snptest
You should have sample files with your phased files, but they have the wrong case-control encoding (they are likely to have plink's 1-2 rather than snptest's 0-1 encoding)

`awk '{if(NR<=2)print; else print $1,$2,$3,$4,$5,$6,$7-1}' study.phased.sample > sample_caco.sample`

NB if you have covariates, don't forget to include those
NB assumes the phenotype is the 7th column and encoded as 1-2; if this is not the case, recode it yourself manually

## 6. run snptest
*6_run_snptest.sh*

NB if your imputed files are compressed (as they should be) you'll need to gzip the sample file you created above as well.

## 7. filter snptest results
*7_filter_snptest.sh*

Filters the snptest output by MAF (case and control separately), HWE (control) and info score. Check the column number are correct (may change depending on your input). Also check the lineSkip field – the new version of snptest adds some header lines.

## 8. summary plots (manhattan + QQ plot)
*8_manhattan_QQ-plots.sh*

I generally just run the R script directly. Also check the beta/se/pvalue columns if they're different.

## 9. meta-analysis
*9_metaLoop.sh*

If you are performing a meta-analysis (with other studies) run meta per chromosome

