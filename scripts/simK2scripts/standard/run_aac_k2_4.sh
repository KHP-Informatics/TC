#!/bin/sh
#$ -N aac_k2_4
#$ -S /bin/sh
#$ -o a_4_out
#$ -e a_4_err
#$ -cwd
#$ -l h_vmem=2G
#$ -V
#$ -t 1-100

Rscript simulate_new_aac.R 4
