#!/usr/bin/perl -w

###########################################
# This program creates a bunch of impute scripts; one for each chunk on each chromosome
# These scripts are have then to be submitted manually to the compute cluster queue
# Heavily edited by Marc Henrion (18 January 2014) from an original script by Sara Dobbins (7 May 2009), and which had already been edited by Yussame Ma (7 July 2010)
# This script now uses both the UK10K (ALSPACE + Twins) and the 100 Genomes data sets as reference panels and merges the reference panel first (option "-merge_ref_panels")
##########################################

use strict;

#-----get input parameters
my $inputPrefix="/scratch/cancgene/plaw/CLL_newdata/phased/CLL_newdata_chr"; # start (incl. of path) of the input files, up to where the chromosome is specified
my $inputSuffix=".phased.haps"; # end of the input filenames, i.e. everything after the chromosome is specified, including the file extension
my $outputDir="/scratch/cancgene/plaw/CLL_newdata/imputed/"; # path to the directory where the output should go
my $outputPrefix="CLL_new_imputed_chr"; # start (excl. of path) of the output files, up to where the chromosome is specified
my $scriptDir="/scratch/cancgene/plaw/CLL_newdata/scripts/"; # directory where the generated script should be written to
my $refPanel1Prefix="/scratch/Dalek/cancgene/studies/imputation_ref_panels/RefPanels_1kG-Dec2013_UK10K-Apr2014_consistent_set_of_variants/1kG_Dec2013_integrated_chr"; # path & filename prefix for the first reference panel files (up to where the chr is specified)
my $refPanel1SuffixHaps=".matched_with_UK10K.hap"; # suffix for the first reference panel haplotype files, including the file extension
my $refPanel1SuffixLeg=".matched_with_UK10K.legend"; # suffix for the first reference panel legend files, including the file extension
my $refPanel2Prefix="/scratch/Dalek/cancgene/studies/imputation_ref_panels/RefPanels_1kG-Dec2013_UK10K-Apr2014_consistent_set_of_variants/UK10K_ALSPAC_TWINS_combined_chr"; # path & filename prefix for the second reference panel files (up to where the chr is specified)
my $refPanel2SuffixHaps=".matched_with_1kG.hap"; # suffix for the second reference panel haplotype files, including the file extension
my $refPanel2SuffixLeg=".matched_with_1kG.legend"; # suffix for the second reference panel legend files, including the file extension
my $geneticMapPrefix="/scratch/Dalek/cancgene/mhenrion/1kG_Dec2013_integrated_unzipped/genetic_map_chr"; # path and filename of the genetic map files, up to where the chromosome is specified
my $geneticMapSuffix="_combined_b37.txt"; # suffix, incl. file extension, of the genetic map files
my @chrSizes=(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566); # the length in bp of the chromosomes that will be run; this is hg19
my $scriptHeader="#!/bin/bash\n#BSUB -n 2\n#BSUB -P prepay-houlston\n#BSUB -W 168:00\n#BSUB -R \"span[hosts=1]\"\n#BSUB -o ".$scriptDir."temp/impute_%J_stdout.txt\n#BSUB -e ".$scriptDir."temp/impute_%J_stderr.txt"; # the first few lines of the submission scripts with queue controller specific settings
my $maxJobsPerNode=2; # maximum number of jobs to be run on any one node
my $chunkSize=5000000; # size of the chunks that will be run in bp
my $imputeProg="/scratch/cancgene/plaw/programs/impute2/impute_v2.3.2_x86_64_static/impute2"; # full path of the imputev2 imputation program

#----- start log file
my $command="mkdir";
system $command, $outputDir unless(-d $outputDir);

my $logfile=$outputDir.$outputPrefix.".log";

open FILELOG, ">", $logfile or die "Cannot open $logfile ($!).";
print FILELOG "This is whole_genome_imputation_RCC_UK10K_1kG_UKGWAS.pl.\n\n";

#-----print parameters to logfile
print FILELOG "Input parameters\n\t";
print FILELOG "inputPrefix = < $inputPrefix >,\n\t";
print FILELOG "inputSuffix = < $inputSuffix >,\n\t";
print FILELOG "refPanel1Prefix = < $refPanel1Prefix >,\n\t";
print FILELOG "refPanel1SuffixHaps = < $refPanel1SuffixHaps >,\n\t";
print FILELOG "refPanel1SuffixLeg = < $refPanel1SuffixLeg >,\n\t";
print FILELOG "refPanel2Prefix = < $refPanel2Prefix >,\n\t";
print FILELOG "refPanel2SuffixHaps = < $refPanel2SuffixHaps >,\n\t";
print FILELOG "refPanel2SuffixLeg = < $refPanel2SuffixLeg >,\n\t";
print FILELOG "geneticMapPrefix = < $geneticMapPrefix >,\n\t";
print FILELOG "geneticMapSuffix = < $geneticMapSuffix >,\n\t";
print FILELOG "outputDir = < $outputDir >,\n\t";
print FILELOG "outputPrefix = < $outputPrefix >,\n\t";
print FILELOG "scriptDir = < $scriptDir >,\n\t";
print FILELOG "imputeProg = < $imputeProg >,\n\t";
print FILELOG "maxJobsPerNode = < $maxJobsPerNode >,\n\t";
print FILELOG "chunkSize = < $chunkSize >,\n\t";
print FILELOG "scriptHeader = < $scriptHeader >,\n\t";
print FILELOG "chrSizes = < @chrSizes >.\n\n";

#-----start loop on chr
for my $chr (1..22){
    #--print info to logfile
    print FILELOG "############\n";
    if($chr>9){
	print FILELOG "## CHR $chr ##\n";
    }else{
	print FILELOG "## CHR 0$chr ##\n";
    }
    print FILELOG "############\n\n";

    #--set start and end positions
    my $usrStart=1;
    my $usrEnd=$chrSizes[$chr-1] + $chunkSize;

    #--input files
    my $geno=$inputPrefix.$chr.$inputSuffix;
    
    my $hapRef1=$refPanel1Prefix.$chr.$refPanel1SuffixHaps;
    my $legRef1=$refPanel1Prefix.$chr.$refPanel1SuffixLeg;
    my $hapRef2=$refPanel2Prefix.$chr.$refPanel2SuffixHaps;
    my $legRef2=$refPanel2Prefix.$chr.$refPanel2SuffixLeg;
    my $map=$geneticMapPrefix.$chr.$geneticMapSuffix;

    #--output files
    my $outputDirTmp=$outputDir."chr".$chr."/";
    my $command="mkdir";
    system $command, $scriptDir unless(-d $scriptDir);
    system $command, $outputDirTmp unless(-d $outputDirTmp);
    system $command, $outputDirTmp."Main/" unless(-d $outputDirTmp."Main");
    system $command, $outputDirTmp."Info/" unless(-d $outputDirTmp."Info");
    system $command, $outputDirTmp."Summary/" unless(-d $outputDirTmp."Summary");
    
    #--generate the commands for the current chromosome
    print FILELOG "generating the commands for chromosome $chr\n";
    
    my @commands=();
    
    my $beginInterval=$usrStart;
    my $endInterval=$usrStart+$chunkSize;
    
    my $countLoops=0;

    while($endInterval<$usrEnd){	
	print FILELOG "generating the commands for the interval from $beginInterval to $endInterval\n";
	
	my $outputMain=$outputDirTmp."Main/".$outputPrefix.$chr."_".$beginInterval."_".$endInterval.".txt";
	my $outputInfo=$outputDirTmp."Info/".$outputPrefix.$chr."_".$beginInterval."_".$endInterval.".txt";
	my $outputSummary=$outputDirTmp."Summary/".$outputPrefix.$chr."_".$beginInterval."_".$endInterval.".txt";

	my $command=$imputeProg." -merge_ref_panels -use_prephased_g -align_by_maf_g -m ".$map." -h ".$hapRef1." ".$hapRef2." -l ".$legRef1." ".$legRef2." -known_haps_g ".$geno." -Ne 20000 -int ".$beginInterval." ".$endInterval." -buffer 250 -o ".$outputMain." -i ".$outputInfo." -r ".$outputSummary;
	
	push @commands, $command;
	
	$countLoops++;
	$beginInterval=$beginInterval+$chunkSize;
	$endInterval=$beginInterval+$chunkSize-1;
	
    }
    print FILELOG "$countLoops loops to be performed for chromosome $chr\n";
    
    my $commandNumber = $#commands+1;
    my $scriptNumber = $commandNumber/$maxJobsPerNode;
    my $num=int($scriptNumber + 0.9999);
    
    #--writing the scripts
    my $n=0;
    for my $i (1..$num){
	my $submitScript=$scriptDir."impute_submit_".$outputPrefix.$chr."_".$i.".sh";
	open FILEOUT, ">", $submitScript or die "Cannot open $submitScript ($!).";
	
	print FILEOUT "$scriptHeader\n";
	
	my $start=$n;
	my $end=$n+$maxJobsPerNode;
	
	for(my $j=$start; $j<$end;$j++){
	    if ($commands[$j]){
		print FILEOUT "$commands[$j] &\n";
	    }
	}
	
	$n=$end;
	
	print FILEOUT "\nwait\n";
	close FILEOUT;

	#--give execution permission to each script
	my $command="chmod";
	push my @args, "u=rwx,g=rx,o=";
	push @args, $submitScript;
	system $command, @args;
    }
    
#-----close the loop over chromosomes
    print FILELOG "Done for chromosome $chr.\n\n";
}

#-----close the log file
print FILELOG "This is the end.\n";
close FILELOG;

