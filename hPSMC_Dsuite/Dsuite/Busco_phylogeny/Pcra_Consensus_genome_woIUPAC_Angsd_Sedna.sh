#!/bin/bash

# this is the header for SEDNA
#SBATCH --job-name=consensus_genome   ## job name
#SBATCH -e consensus_genome%j.e.txt    ## error message name
#SBATCH -o consensus_genome.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -p standard
#SBATCH -c 5    ## <number of cores to ask for> (cpu cores per task?)
#SBATCH --array=1-7 # -81%30 
#SBATCH --mem=20G # max memory on standard nodes = 96G; himem = 1.5TB
#SBATCH -t 96:00:00   ## walltime in mins or mins:secs or hrs:mins:secs. 
#  SBATCH -D /scratch/pmorin/temp		## <folder to change to when starting the job>

##########################################################################################
# Uses Angsd to extract consensus genome

# run from scripts folder (scripts/03_ANGSD)
module load bio/angsd/0.933
# module load bio/samtools/1.11 # to index bam files prior to running ANGSD
# module load tools/rclone

set -eux

BAMLIST=(
Morgan_Oorc_dedup.bam
z0018462_Pcra_dedup.bam
z0027510_Pcra_dedup.bam
z0045928_Pcra_dedup.bam
z0112653_Gmac_dedup.bam
z0112731_Gmac_dedup.bam
z0017981_Gmac_dedup.bam
)
# BAMLIST=Delphinid_bam_file_paths_test.txt # set array number to match number of bam files
REFDIR=/home/pmorin/Ref_genomes/Oorc/NCBI_GCA_937001465.1
REF=GCA_937001465.1_mOrcOrc1.1_genomic.fna
THREADS=10

# input and output directories
INDIR=/home/pmorin/projects/Pcra/hPSMC_Dsuite_June2025/Dsuite # input directory
OUTDIR=${INDIR}/Consensus_wo_iupac # output directory
TEMP=${OUTDIR}/TEMP
mkdir -p ${OUTDIR}
mkdir -p ${TEMP}
baseQ=30
mapQ=25

#####################
# Get working bam file based on array number
NUM=$(printf %02d ${SLURM_ARRAY_TASK_ID})

# BAMFILE
BAMFILE=${BAMLIST[$SLURM_ARRAY_TASK_ID-1]} # if using an alias list within the script
echo ${BAMFILE}

SP=`echo $BAMFILE | cut -f2 -d "_"`
sample=`echo $BAMFILE | cut -f1 -d "_"`

# Define input and output file names for ANGSD
REF_GENOME=${REFDIR}/${REF}
OUT_PREFIX=${SP}_${sample}_angsd_consensus_genome

#########
# Create consensus genome
angsd -i ${INDIR}/${BAMFILE} -ref ${REF_GENOME} -out ${OUTDIR}/$OUT_PREFIX -minQ $baseQ -minMapQ $mapQ -doMajorMinor 1 -doMaf 1 -doCounts 1 -GL 1 -doFasta 2 -doDepth 1

# mv *genome.arg ${TEMP}
mv ${OUTDIR}/*genome.mafs.gz ${TEMP}/

# copy consensus sequence to google Drive with Rclone
# rclone copy -P ${BW_CON}/$OUTPUT_PREFIX.fa.gz Gdrive:

# Print job completion time
date

# -doFasta 2   use the most common base. In the case of ties a random base is chosen among the bases with the same maximum counts. N's or filtered based are ignored. The "-doCounts 1" options for allele counts is needed in order to determine the most common base. If multiple individuals are used the four bases are counted across individuals.
