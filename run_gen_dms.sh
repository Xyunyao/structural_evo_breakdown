#!/bin/bash

# Load conda environment
source ~/miniconda3/etc/profile.d/conda.sh   # adjust this path if needed
conda activate struct-evo

# Optionally add repo to PYTHONPATH (if dms_utils isn’t installed in env)
export PYTHONPATH="/mnt/beegfs/home/friesner/yx2218/structural-evolution/bin:$PYTHONPATH"

# Input/output paths
fastafile=/mnt/beegfs/home/friesner/yx2218/structural-evolution/structural_evo_breakdown/SHP2/SHP2_FULL.fasta
output_dir=./SHP2/DMS_library/SHP2_FULL_dms.fasta

# Run script
python ./generate_dms_fasta.py --protein_ID SHP2 --fastafile "$fastafile" --dmspath "$output_dir"

echo "DMS fasta file generated at: $output_dir"

# run_dm_job() {
#     local seq_path=$1
#     local output_path=$2

#     # Run the DM job
#     echo "Running DM job with sequence file: $seq_path"
#       # Simulate some processing time
#     echo "DM job completed. Output written to: $output_path"
# }
