#!/bin/bash -l

#$ -N 		## Name the job
#$ -j y         ## Merge error & output files
#$ -pe omp 4	## Request resources

## This script is based on a university computing cluster and is likely to
## vary based on your computing environment
## Learn more about high performance computing:
##		https://github.com/Climate-CAFE/hpc_batch_jobs_micro_tutorial

module load R/4.4.0
module load libsodium/1.0.18
Rscript /pathtoscript/04A_UTCI_WBGT_Intermediate_Dynamic_v2.R $SGE_TASK_ID

## Run from b = 2 (because we add the previous month)
## In Terminal, cd to the directory in which this bash script is located. 
## qsub -P project -t 2-25 04B_UTCI_WBGT_Intermediate_Dynamic_v2.sh
##
## This script is set to run as an "array". The indices from 2 to 25
## will be passed to the code so that each month of ERA5 data is 
## processed separately. 