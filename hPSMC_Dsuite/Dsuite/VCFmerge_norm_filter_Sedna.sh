#!/bin/bash

# this is the header for SEDNA
#SBATCH --job-name=VCFfilter   ## job name
#SBATCH -e VCFfilter%j.e.txt    ## error message name
#SBATCH -o VCFfilter.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -c 10    ## (cpu cores per task)
#SBATCH -p medmem  #partitians= standard= 96G; medmem=192G; himem = 1.5TB
#  SBATCH --mem=80G # max memory
#SBATCH -t 4-0   ## walltime in mins, mins:secs, hrs:mins:secs, days-hours 
#SBATCH --ntasks=1                   ## Run a single task
#  SBATCH --array=1-14 # 
#  SBATCH -D /scratch/pmorin/temp		## <folder to change to when starting the job>

##########################################################################################

# extract variable sites from bam file to VCF file.

# module load bio/samtools/1.19
module load bio/bcftools/1.11

# INDIR=/home/pmorin/projects/Miscellaneous/Ziphiidae_ABBA_BABA/Mesoplodon_map2Bbai/vcf_files
OUTDIR=/home/pmorin/projects/Pcra/hPSMC_Dsuite_June2025/Dsuite/vcf_files
TEMPDIR=/scratch/pmorin/temp/Pcra_VCF_merge
REF=/home/pmorin/Ref_genomes/Oorc/NCBI_GCA_937001465.1/GCA_937001465.1_mOrcOrc1.1_genomic.fna
THREADS=10

mkdir -p ${OUTDIR}
mkdir -p ${TEMPDIR}

############### 
# merge, normalize, and filter vcf files
bcftools merge -Oz -o ${TEMPDIR}/Ziphiid_joint_variants.vcf.gz ${TEMPDIR}/*vcf.gz

# normalize indels and index:
bcftools norm -m -any -f ${REF} -Oz -o ${TEMPDIR}/Ziphiid_normalized.vcf.gz ${TEMPDIR}/Ziphiid_joint_variants.vcf.gz

bcftools index --threads ${THREADS} ${TEMPDIR}/Ziphiid_normalized.vcf.gz

# filter variants
bcftools filter --threads ${THREADS} -i 'QUAL>30 && FORMAT/DP>10' -Oz ${TEMPDIR}/Ziphiid_normalized.vcf.gz -o ${OUTDIR}/Ziphiid_filtered.vcf.gz 

# bcftools output formats: -O (VCF: z=compressed, v=uncompressed. BCV: b=compressed, u=uncompressed)


