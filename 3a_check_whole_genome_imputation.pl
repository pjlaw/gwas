#!/usr/bin/perl

#-----get input parameters
my $study="CLL_new_imputed_chr";
#my $caco="CO_P";
#my $outputPrefix=$study."_".$caco."_uk10k_1kG_chr"; # start (excl. of path) of the output files from impute, up to where the chromosome is specified
#my $outputDir="/scratch/cancgene/plaw/crc/uk10k_1kG_impute/".$study."/"; # path to the directory where the output from impute went
my $outputDir="/scratch/cancgene/plaw/CLL_newdata/imputed/"; # path to the directory where the output from impute went
my $outputPrefix=$study; # start (excl. of path) of the output files from impute, up to where the chromosome is specified
my $summaryLog=$study."_impute_check_chr1-22.log"; # the name of the log file that will be generated and contains the stats of which chunks have failed / are missing
my @chrSizes=(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566); # the length in bp of the chromosomes that will be run; this is hg19
my $chunkSize=5000000; # size of the chunks that have been run in bp

#-----start log file
my $logfile=$outputDir.$summaryLog;

open FILELOG, ">", $logfile or die "Cannot open $logfile ($!).";
print FILELOG "This is check_whole_genome_imputation_RCC_UK10K_1kG_NCIGWAS.pl.\n\n";
print FILELOG "Input parameters:\n\t";
print FILELOG "ouputDir = < $outputDir >,\n\t";
print FILELOG "outputPrefix = < $outputPrefix >,\n\t";
print FILELOG "chunkSize = < $chunkSize >.\n\n";

#-----loop over chromosomes, checking that output files for each chunk exist and, if not, if this is properly reported in the summary files
for my $chr (1..22){

    print FILELOG "########################################\n";
    print FILELOG "## Checking files for chromosome $chr. ##\n";
    print FILELOG "########################################\n\n";

    #--loop over the chunks
    my $usrStart=1;
    my $usrEnd=$chrSizes[$chr-1] + $chunkSize;
    
    my $beginInterval=$usrStart;
    my $endInterval=$usrStart+$chunkSize;

    while($endInterval<$usrEnd){	
	if(-e $outputDir."chr$chr/Main/$outputPrefix".$chr."_".$beginInterval."_".$endInterval.".txt"){
	    my $sumFile=$outputDir."chr$chr/Summary/$outputPrefix".$chr."_".$beginInterval."_".$endInterval.".txt";		
	    my $accCheck=qx/tail -1 $sumFile/;
	    $accCheck=~s/^\s+//;
	    my @accCheck=split /\s+/, $accCheck;
	    my $accChunk=$accCheck[2];
	    if ($accChunk eq "were" || $accChunk eq "was"){
			my $numWarnings = $accCheck[0];
			my $accCheck=qx/tail -3 $sumFile | head -1 /;
			$accCheck=~s/^\s+//;
			my @accCheck=split /\s+/, $accCheck;
			my $accChunk=$accCheck[2];
			if($accChunk>95){
				print FILELOG "SUPER2: output for chr $chr, interval $beginInterval -> $endInterval exists and imputation accuracy is < $accChunk > BUT $numWarnings warnings.\n";
			}else{
				print FILELOG "CAUTION2: output for chr $chr, interval $beginInterval -> $endInterval exists, but imputation accuracy is only < $accChunk > with $numWarnings warnings.\n";
			}
		} else{
	    if($accChunk>95){
			print FILELOG "SUPER: output for chr $chr, interval $beginInterval -> $endInterval exists and imputation accuracy is < $accChunk >.\n";
		}else{
			print FILELOG "CAUTION: output for chr $chr, interval $beginInterval -> $endInterval exists, but imputation accuracy is only < $accChunk >.\n";
			}
		}
	}else{
	    my $sumFile=$outputDir."chr$chr/Summary/$outputPrefix".$chr."_".$beginInterval."_".$endInterval.".txt";		
	    if(-e $sumFile){
			my $commandOutput=qx/tail -2 $sumFile | head -1 | grep "There are no SNPs in the imputation interval" | wc -l/;
			chomp $commandOutput;
			if($commandOutput eq "1"){
				print FILELOG "OK: chr $chr, interval $beginInterval -> $endInterval had no SNPs in the imputation interval and reported it so.\n";
			}else{
				my $commandOutput=qx/tail -2 $sumFile | head -1 | grep "there will not be any SNPs in the output file" | wc -l/;
				chomp $commandOutput;
				if($commandOutput eq "1"){
				print FILELOG "OK: chr $chr, interval $beginInterval -> $endInterval had no SNPs in the imputation interval and reported it so.\n";
				}else{
				print FILELOG "WARNING: chr $chr, interval $beginInterval -> $endInterval seems to have finished without leaving a reason why. Summary file < $sumFile > exists but possibly left unfinished (?).\n";
				}
		}
	    }else{
			print FILELOG "WARNING: chr $chr, interval $beginInterval -> $endInterval seems to have finished without leaving a reason why (no summary file at < $sumFile >).\n";
	    }
	}
	
	$beginInterval=$beginInterval+$chunkSize;
	$endInterval=$beginInterval+$chunkSize-1

    #--close the loop over chunks
    }

    print FILELOG "Done.\n\n"

#-----close the loop over chromomes
}

#-----exit
print FILELOG "This is the end.\n";
close FILELOG;
