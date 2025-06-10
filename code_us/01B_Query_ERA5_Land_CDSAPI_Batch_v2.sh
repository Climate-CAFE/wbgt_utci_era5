#!/bin/bash -l

#$ -N 			## Name the job
#$ -j y         	## Merge error & output files
#$ -pe omp 16		## Request resources
#$ -l h_rt=36:00:00	## Set runtime

## This script is based on a university computing cluster and is likely to
## vary based on your computing environment
## Learn more about high performance computing:
##		https://github.com/Climate-CAFE/hpc_batch_jobs_micro_tutorial

module load R/4.4.0
module load libsodium/1.0.18
Rscript /pathtoscript/01A_Query_ERA5_Land_CDSAPI_Batch_v2.R

## In Terminal, cd to the directory in which this bash script is located. 
## qsub -P project 01B_Query_ERA5_Land_CDSAPI_Batch_v2.sh
