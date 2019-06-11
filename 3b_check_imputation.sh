#!/bin/bash

path=/scratch/cancgene/plaw/CLL2_X/imputed
fname=$path/CLL1_impute_check_chrX.log

grep SUPER $fname > $path/ok.log
grep OK $fname >> $path/ok.log
grep CAUTION $fname > $path/caution.log
grep WARNING $fname > $path/warning.log
