#!/bin/bash
#submit with  ./submit_workflow.sh <inputDir>
set -e

INPUT_DIR=$1

SCRIPTS_DIR=/grphome/grp_Assembly/pcrScripts/snakemake/scripts

# Create a unique tmux session name based on the input directory
SESSION_NAME="snakemake_$(basename "$INPUT_DIR")_$(date +%Y%m%d_%H%M%S)"
LOGFILE="snakemake_$(date +%Y%m%d_%H%M%S).log"

# Load modules and run snakemake in a detached tmux session
tmux new-session -d -s "$SESSION_NAME" bash -c "
    module load miniconda3
    conda run -p /grphome/grp_Assembly/snakemake --no-capture-output \
        snakemake --executor slurm --jobs 200 \
            --latency-wait 60 \
            --use-conda \
            --default-resources mem_mb_per_cpu=4000 \
            -s main.smk \
            --config input_dir='$INPUT_DIR' scripts_dir='$SCRIPTS_DIR' \
            2>&1 | tee $LOGFILE
    echo 'Workflow completed. Press enter to close this tmux session.'
    read
"

echo "============================================"
echo "Snakemake workflow launched in tmux session: $SESSION_NAME"
echo "Log file: $LOGFILE"
echo ""
echo "To attach and view progress:"
echo "  tmux attach -t $SESSION_NAME"
echo ""
echo "To detach (once attached):"
echo "  Press Ctrl+b, then d"
echo ""
echo "To kill the session:"
echo "  tmux kill-session -t $SESSION_NAME"
echo "============================================"
