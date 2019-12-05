#!/bin/bash
#BSUB -J plots
#BSUB -n 1
#BSUB -W 120:00
#BSUB -P prepay-houlston
#BSUB -o plots_%J_stdout.txt
#BSUB -e plots_%J_stderr.txt

R --vanilla < 8_manhattan_QQ-plots.R 

