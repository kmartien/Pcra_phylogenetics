#!/bin/bash

# this is the header for SEDNA
#SBATCH --job-name=bam2vcf   ## job name
#SBATCH -e bam2vcf%j.e.txt    ## error message name
#SBATCH -o bam2vcf.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -c 5    ## (cpu cores per task)
#SBATCH -p medmem  #partitians= standard= 96G; medmem=192G; himem = 1.5TB
#  SBATCH --mem=40G # max memory
#SBATCH -t 4-0   ## walltime in mins, mins:secs, hrs:mins:secs, days-hours 
#SBATCH --ntasks=1                   ## Run a single task
#SBATCH --array=1-7 # 
# SBATCH -D /scratch/pmorin/temp		## <folder to change to when starting the job>

##########################################################################################

# extract variable sites from bam file to VCF file.

# module load bio/samtools/1.19
module load bio/bcftools/1.11

INDIR=/home/pmorin/projects/Pcra/hPSMC_Dsuite_June2025/Dsuite
OUTDIR=/home/pmorin/projects/Pcra/hPSMC_Dsuite_June2025/Dsuite/vcf_files
TEMPDIR=/scratch/pmorin/temp/Pcra_VCF_merge
REF=/home/pmorin/Ref_genomes/Oorc/NCBI_GCA_937001465.1/GCA_937001465.1_mOrcOrc1.1_genomic.fna
BAMLIST=${INDIR}/Pcra_VCF_bamlist.txt

mkdir -p ${OUTDIR}
mkdir -p ${TEMPDIR}


# Get working bam file based on array number
NUM=$(printf %02d ${SLURM_ARRAY_TASK_ID})

# generate windows list for each scaffold
BAM=$(head -n ${NUM} ${BAMLIST} | tail -n 1)
echo ${BAM}

filename=${BAM} # `echo $BAM | cut -f7 -d "/"`
sampleID=`echo $filename | cut -f1 -d "_"`
species=`echo $BAM | cut -f2 -d "_"`

bcftools mpileup -B -q 30 -Q 30 -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR -f ${REF} -Ou --threads 20 ${BAM} \
 | bcftools call -Oz --threads 20 --multiallelic-caller -o ${TEMPDIR}/${species}_${sampleID}.vcf.gz
# -a flag is for various quality score headers in the BCF output. I don't use them at the moment but it was something standard that I had seen others add, so I included those -a flags in case it's needed for later (Morgan)
# -Oz = compressed vcf format (vcf.gz); -Ob=compressed bcf format (bcf.gz)

# bcftools view "${species}_${sampleID}.bcf" -o "${species}_${sampleID}.vcf"

bcftools index --threads 20 ${TEMPDIR}/${species}_${sampleID}.vcf.gz




