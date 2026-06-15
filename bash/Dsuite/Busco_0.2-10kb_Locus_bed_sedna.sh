#!/bin/bash
#SBATCH --job-name=Busco_bed     ## job name
#SBATCH -e Busco_bed_%j.e.txt    ## error message name
#SBATCH -o Busco_bed.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -p medmem
#SBATCH -c 1
#SBATCH --array=1-35 # -81%20 # %# limits the number of subjobs in an array that run concurrently, e.g., array 1-10, with 2 jobs at a time, would be --array=1-10%2.
#  SBATCH --mem=4G
#SBATCH -t 2-0   ## walltime in mins or mins:secs or hrs:mins:secs. 
# SBATCH -D /scratch/pmorin/temp		## <folder to change to when starting the job>

######################################################################
# use output from BUSCO to generate a bed file of 'complete' genes with size range from 200bp to 10kb.

INDIR=
OUTDIR=
INFILE=
Min=200
Max=15000
