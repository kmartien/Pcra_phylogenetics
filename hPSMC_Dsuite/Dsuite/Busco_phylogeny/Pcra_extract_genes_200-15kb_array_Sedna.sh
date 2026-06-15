#!/bin/bash
#SBATCH --job-name=extract_seqs     ## job name
#SBATCH -e extract_seqs_%j.e.txt    ## error message name
#SBATCH -o extract_seqs.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -p medmem
#SBATCH -c 1
#SBATCH --array=1-35 # -81%20 # %# limits the number of subjobs in an array that run concurrently, e.g., array 1-10, with 2 jobs at a time, would be --array=1-10%2.
#  SBATCH --mem=4G
#SBATCH -t 2-0   ## walltime in mins or mins:secs or hrs:mins:secs. 
# SBATCH -D /scratch/pmorin/temp		## <folder to change to when starting the job>

######################################################################

# Load necessary modules
module load bio/bedtools/2.31.1
module load bio/samtools/1.19
module load bio/htslib/1.19
# module load tools/pigz/2.4

# Define input and output directories and files
THREADS=5
PROJDIR=/home/pmorin/projects/Pcra/hPSMC_Dsuite_June2025/Dsuite
FASTADIR=$PROJDIR/Consensus_wo_iupac
# REFDIR=/home/pmorin/Ref_genomes/Oorc/NCBI_GCA_937001465.1

BEDFILE=$PROJDIR/Oorc_BCA_937001465.1_Busco_200-15k.bed
OUTPUTDIR=$PROJDIR/"Busco_loci_out_200-15k"
mkdir -p $OUTPUTDIR

FASTAFILES=(
Gmac_z0017981_angsd_consensus_genome.fa.gz
Gmac_z0112653_angsd_consensus_genome.fa.gz
Gmac_z0112731_angsd_consensus_genome.fa.gz
Oorc_Morgan_angsd_consensus_genome.fa.gz
Pcra_z0018462_angsd_consensus_genome.fa.gz
Pcra_z0027510_angsd_consensus_genome.fa.gz
Pcra_z0045928_angsd_consensus_genome.fa.gz
)

cd $FASTADIR

# Extract sequence based on array index
FASTA=${FASTAFILES[$SLURM_ARRAY_TASK_ID-1]}

FA=`echo $FASTA | cut -f1 -d "."` # fasta file name before ".fa.gz"
SP=`echo $FASTA | cut -f1 -d "_"`
SPID=`echo $FASTA | cut -f1,2 -d "_"`

echo "indexing $FASTA..."
# 
# # unzip fasta files # this shouldn't be necessary with bedtools 2.31.1 (fixed issue)
# # gunzip $FASTADIR/$FASTA
# 
# index fasta files
samtools faidx $FASTADIR/$FASTA

# echo "Extracting sequences from $FASTA..."

# Use bedtools to extract sequences from fasta file based on bed file coordinates
bedtools getfasta -fi $FASTADIR/$FASTA -bed $BEDFILE -fo $OUTPUTDIR/${FA%}.bedout2.fa -name+
# -name+ adds information from additional columns to the fasta headers

# add species name to each fasta header
sed -E "s/^>/>${SPID}_/g" $OUTPUTDIR/${FA%}.bedout2.fa > $OUTPUTDIR/${FA%.bam}.bedout.fa

# rm $OUTPUTDIR/${FA%.fa}.bedout2.fa # for some reason this causes some to fail, so leave off of script and remove manually

