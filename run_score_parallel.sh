#!/bin/bash
#SBATCH --job-name=dms_score
#SBATCH --partition=fcpu        # CPU partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=0-1              # adjust to number of lines in pdb_jobs.txt
#SBATCH --output=logs/dms_%A_%a.out
#SBATCH --error=logs/dms_%A_%a.err

set -euo pipefail

# Load conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate struct-evo

# Add repo to PYTHONPATH
export PYTHONPATH="/mnt/beegfs/home/friesner/yx2218/structural-evolution/bin:${PYTHONPATH:-}"

# Paths
fastafile="/mnt/beegfs/home/friesner/yx2218/structural-evolution/structural_evo_breakdown/SHP2/SHP2_FULL.fasta"
dmslib_dir="./SHP2/DMS_library"
output_dir="./SHP2/output"
mkdir -p "$dmslib_dir" "$output_dir" logs

# Pick job line based on SLURM_ARRAY_TASK_ID
jobline=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" pdb_jobs2.txt)
pdbfile=$(echo "$jobline" | awk '{print $1}')
chain=$(echo "$jobline" | awk '{print $2}')
start_index=$(echo "$jobline" | awk '{print $3}')

# Output file with timestamp
outfile="$output_dir/$(basename "$pdbfile" .pdb)_chain${chain}_start${start_index}_dms_scores_$(date +%Y%m%d_%H%M%S).csv"

echo "[$(date)] Running job $SLURM_ARRAY_TASK_ID on $pdbfile (chain $chain, start_index $start_index)..."

python /mnt/beegfs/home/friesner/yx2218/structural-evolution/structural_evo_breakdown/score_log_likelihoods_v3.py \
    --pdbfile "$pdbfile" \
    --start_index "$start_index" \
    --seqfile_path "$fastafile" \
    --seqpath "$dmslib_dir/SHP2_FULL_dms.fasta" \
    --chain "$chain" \
    --outpath "$outfile" \
    --singlechain-backbone

echo "[$(date)] Job $SLURM_ARRAY_TASK_ID finished. Output: $outfile"
