#!/bin/bash

# Load conda environment
source ~/miniconda3/etc/profile.d/conda.sh   # adjust this path if needed
conda activate struct-evo

# Optionally add repo to PYTHONPATH (if dms_utils isn’t installed in env)
export PYTHONPATH="/mnt/beegfs/home/friesner/yx2218/structural-evolution/bin:$PYTHONPATH"

# Input/output paths
fastafile=/mnt/beegfs/home/friesner/yx2218/structural-evolution/structural_evo_breakdown/SHP2/SHP2_FULL.fasta
output_dir=./SHP2/DMS_library/SHP2_FULL_dms.fasta

# Run script with a function with options
run_dm_job() {

    python ./score_log_likelihood_v3.py \
            --pdbfile $1 \
            --start_index $2  \
            --seqfile_path $3 \
            --seqpath $4 \
            --chain $5 \
            --outpath $6 \
            --singlechain-backbone
    # Simulate some processing time
    echo "DM job completed."
}

run_dm_job 
    "/mnt/beegfs/home/friesner/yx2218/structural-evolution/structural_evo_breakdown/SHP2/2SHP.pdb" \
    1 \
    $fastafile \
    $output_dir \
    "A" \
    ""

# G6 Light Chain
run_protein_mpnn \
    "data/ab_mutagenesis_expts/g6/g6_2fjg_lc_lib.fasta" \
    "L" \
    "output/ab_mutagenesis_expts/g6/proteinMpnnLC" \
    "data/ab_mutagenesis_expts/g6/2fjg_vlh_fvar.pdb"


echo "DMS AI score finished"




