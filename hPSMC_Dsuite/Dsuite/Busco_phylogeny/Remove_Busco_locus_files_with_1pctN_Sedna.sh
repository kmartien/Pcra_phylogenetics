#!/bin/bash
#SBATCH --job-name=remove_pctN
#SBATCH -e remove_pctN_%j.e.txt    ## error message name
#SBATCH -o remove_pctN.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -p medmem
#SBATCH -c 1
# SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --output=fasta-processing.%j.out

# Set path to directory containing the fasta files
FASTA_DIR=/home/pmorin/projects/Pcra/hPSMC_Dsuite_June2025/Dsuite/aligned_loci_out

# Set path to output directory
OUTPUT_DIR=${FASTA_DIR}/loci_high_Ns

mkdir -p ${OUTPUT_DIR}

pctN=0.1

# Change to the fasta directory
cd $FASTA_DIR

# Loop through all fasta files in the directory
for FILE in *.fasta
do
    # Get the total length of all sequences in the file
    LENGTH=$(grep -v '>' $FILE | tr -d '\n' | wc -c)
    
    # Get the number of N's in the file
    NUM_NS=$(grep -v '>' $FILE | tr -d '\n' | tr -cd 'Nn' | wc -c)
    
    # Calculate the percentage of N's in the file
    PERCENT_NS=$(echo "scale=4; ($NUM_NS/$LENGTH)*100" | bc)
    
    # If the percentage of N's is greater than pctN, rename the file
    if (( $(echo "$PERCENT_NS > $pctN" | bc -l) )); then
        mv $FILE ${OUTPUT_DIR}/${FILE}_high_N.fasta
    fi
done

# move loci with <1% N's to new folder
lowN=Loci_0.1pct_low_Ns
mkdir ${FASTA_DIR}/${lowN}
mv ${FASTA_DIR}/*fasta ${lowN}
