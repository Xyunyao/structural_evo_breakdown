#!/bin/bash
set -euo pipefail

# Load conda environment
source ~/miniconda3/etc/profile.d/conda.sh   # adjust path if needed
conda activate struct-evo

# Optionally add repo to PYTHONPATH (if dms_utils isn’t installed in env)
export PYTHONPATH="/mnt/beegfs/home/friesner/yx2218/structural-evolution/bin:${PYTHONPATH:-}"


# Input/output paths
fastafile="/mnt/beegfs/home/friesner/yx2218/structural-evolution/structural_evo_breakdown/SHP2/SHP2_FULL.fasta"
dmslib_dir="./SHP2/DMS_library"
mkdir -p "$dmslib_dir"

# Run script with a function with options
run_dm_job() {
    local pdbfile="$1"
    local start_index="$2"
    local seqfile_path="$3"
    local seqpath="$4"
    local chain="$5"
    local outpath="$6"

    echo "Running DM job on $pdbfile (chain $chain)..."
    python /mnt/beegfs/home/friesner/yx2218/structural-evolution/structural_evo_breakdown/score_log_likelihoods_v3.py \
        --pdbfile "$pdbfile" \
        --start_index "$start_index" \
        --seqfile_path "$seqfile_path" \
        --seqpath "$seqpath" \
        --chain "$chain" \
        --outpath "$outpath" \
        --singlechain-backbone

    echo "DM job completed for $pdbfile."
}

# Test run
run_dm_job \
    "/mnt/beegfs/home/friesner/yx2218/structural-evolution/structural_evo_breakdown/SHP2/structure_data/4DGP.pdb" \
    0 \
    "$fastafile" \
    "$dmslib_dir/SHP2_FULL_dms.fasta" \
    "A" \
    "./SHP2/output/4DGP_dms_scores_$(date +%Y%m%d_%H%M%S).csv"

echo "DMS AI score finished ✅"





