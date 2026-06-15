#!/bin/bash

# this is the header for SEDNA
#SBATCH --job-name=VCFsplit_chrom   ## job name
#SBATCH -e VCFsplit_chrom_%j.e.txt    ## error message name
#SBATCH -o VCFsplit_chrom.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -c 20    ## (cpu cores per task)
#SBATCH -p standard  #partitians= standard= 96G; medmem=192G; himem = 1.5TB
#SBATCH --mem=80G # max memory
#SBATCH -t 3-0   ## walltime in mins, mins:secs, hrs:mins:secs, days-hours 
#SBATCH --ntasks=1                   ## Run a single task
#  SBATCH --array=1 # 
#  SBATCH -D /scratch/pmorin/temp		## <folder to change to when starting the job>

##########################################################################################
# module load GCC/13.2.0 BCFtools/1.19
module load bio/bcftools/1.11

# NOTE: Ensure your VCF file with all chromosomes is BGZF (bgzip) compressed!

VCFDIR=/home/pmorin/projects/Pcra/hPSMC_Dsuite_June2025/Dsuite/vcf_files
VCF=${VCFDIR}/Pcra_Ziphiid_filtered.vcf.gz
OUT="split-vcf/"

# Create new output directory
mkdir -p $OUT

# Ensure the VCF is bgzipped and indexed
# if [ ! -f "${VCF}.tbi" ] && [ ! -f "${VCF}.csi" ]; then
#     echo "Indexing $VCF..."
#     bcftools index "$VCF"
# fi

# Get contig/chromosome names from the VCF header
chroms=$(bcftools view -h "$VCF" | grep '^##contig' | sed -E 's/.*ID=([^,>]+).*/\1/')

# Split the VCF
for chrom in $chroms; do
    echo "Extracting $chrom..."
    bcftools view --threads 20 -r "$chrom" "$VCF" -Oz -o "${OUT}/${chrom}.vcf.gz"
    bcftools index --threads 20 "${OUT}/${chrom}.vcf.gz"
done
