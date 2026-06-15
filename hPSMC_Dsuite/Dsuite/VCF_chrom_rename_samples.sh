#!/bin/bash

# this is the header for SEDNA
#SBATCH --job-name=VCF_rename   ## job name
#SBATCH -e VCF_rename_%j.e.txt    ## error message name
#SBATCH -o VCF_rename.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -c 2    ## (cpu cores per task)
#SBATCH -p standard  #partitians= standard= 96G; medmem=192G; himem = 1.5TB
#SBATCH --mem=10G # max memory
#SBATCH -t 1-0   ## walltime in mins, mins:secs, hrs:mins:secs, days-hours 
#SBATCH --ntasks=1                   ## Run a single task
#  SBATCH --array=1-22 # 
#  SBATCH -D /scratch/pmorin/temp		## <folder to change to when starting the job>

##########################################################################################
# change species/sample names in merged vcf file

module load bio/bcftools/1.11

INDIR=/home/pmorin/projects/Miscellaneous/Ziphiidae_ABBA_BABA/Ziph_outgroup_Dsuite/vcf_files/split-vcf
OUTDIR=${INDIR}
mkdir -p ${OUTDIR}

Rename_file=Hpla_rename_sample.txt

cd ${INDIR}

for file in *vcf.gz
do

shortID=`echo $file | cut -f1,2 -d "."`

# rename samples in vcf files
bcftools reheader -s ${INDIR}/${Rename_file} ${INDIR}/${file} > ${OUTDIR}/${shortID}_renamed.vcf.gz

# need to index the files after renaming
bcftools index ${OUTDIR}/${shortID}_renamed.vcf.gz

done