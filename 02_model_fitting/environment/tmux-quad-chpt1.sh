#!/usr/bin/env bash
# tmux-quad-panes.sh
# Creates a 2x2 tmux layout for R analysis workflow
# Pane 0 (top-left): htop
# Pane 1 (top-right): R console with all modules/env vars loaded
# Pane 2 (bottom-left): GPU monitor
# Pane 3 (bottom-right): working shell
#
# Usage:
#   ./tmux-quad-panes.sh [SESSION_NAME]
#   WORKDIR=~/proj ./tmux-quad-panes.sh chpt1

set -euo pipefail

SESSION="${1:-chpt1}"

# ---- User-tunable via env vars (defaults) ----
WORKDIR="${WORKDIR:-/gpfs/data1/vclgp/lmaden/chpt1}"
CONDA_ENV="${CONDA_ENV:-}"
GPU_WATCH="${GPU_WATCH:-auto}"

# Expand leading ~ in WORKDIR if present
case "$WORKDIR" in
  "~" | "~/"* ) WORKDIR="${WORKDIR/#\~/$HOME}";;
esac

have() { command -v "$1" >/dev/null 2>&1; }

attach_or_switch() {
  if [ -n "${TMUX:-}" ]; then
    exec tmux switch-client -t "$SESSION"
  else
    exec tmux attach -t "$SESSION"
  fi
}

if tmux has-session -t "$SESSION" 2>/dev/null; then
  # Session exists: jump to it
  attach_or_switch
fi

# Create session detached, window 0 named "quad"
tmux new-session -d -s "$SESSION" -n quad -c "$WORKDIR"

# Build a 2x2 grid
tmux split-window -h -t "$SESSION":0
tmux select-pane -L -t "$SESSION":0
tmux split-window -v -t "$SESSION":0
tmux select-pane -R -t "$SESSION":0
tmux split-window -v -t "$SESSION":0
tmux select-layout -t "$SESSION":0 tiled

# Quality-of-life settings for this session only
tmux set-option -t "$SESSION" mouse on
tmux set-option -t "$SESSION" history-limit 50000

# Pane IDs (typical numbering order)
P0="$SESSION:0.0"   # top-left
P1="$SESSION:0.1"   # top-right
P2="$SESSION:0.2"   # bottom-left
P3="$SESSION:0.3"   # bottom-right

# ---- Pane 0: htop/top ----
if have htop; then
  tmux send-keys -t "$P0" "htop" C-m
else
  tmux send-keys -t "$P0" "top" C-m
fi

# ---- Pane 1: R Console with full environment ----
# First, set up environment
tmux send-keys -t "$P1" "cd \"$WORKDIR\"" C-m
tmux send-keys -t "$P1" "module purge" C-m
tmux send-keys -t "$P1" "module load rh9/R/4.5.0" C-m
tmux send-keys -t "$P1" "module load rh9/gdal/3.11.0" C-m
tmux send-keys -t "$P1" "module load rh9/geos 2>/dev/null || true" C-m
tmux send-keys -t "$P1" "module load rh9/proj 2>/dev/null || true" C-m
tmux send-keys -t "$P1" "module load rh9/udunits 2>/dev/null || true" C-m

# Set environment variables
tmux send-keys -t "$P1" "export PROJECT_ROOT=/gpfs/data1/vclgp/lmaden/chpt1" C-m
tmux send-keys -t "$P1" "export DATA_CSV='sampled_data_FINAL_with_lc.csv'" C-m
tmux send-keys -t "$P1" "export READ_STRATEGY=chunked_sample" C-m
tmux send-keys -t "$P1" "export CHUNK_ROWS=400000" C-m
tmux send-keys -t "$P1" "export CHUNK_SAMPLE_FRAC=0.05" C-m
tmux send-keys -t "$P1" "export EDA_ENABLE=TRUE" C-m
tmux send-keys -t "$P1" "export EDA_SAMPLE_FRAC=0.15" C-m
tmux send-keys -t "$P1" "export EDA_MAX_N=300000" C-m
tmux send-keys -t "$P1" "export EDA_MAX_PER_GROUP=40000" C-m
tmux send-keys -t "$P1" "export CPU_FRAC=0.09" C-m
tmux send-keys -t "$P1" "export RAM_FRAC=0.09" C-m
tmux send-keys -t "$P1" "export TMPDIR=/gpfs/data1/vclgp/lmaden/chpt1/tmp" C-m
tmux send-keys -t "$P1" "export IMG_META_WORKERS=30" C-m

# Print welcome message
tmux send-keys -t "$P1" "echo '==========================================='" C-m
tmux send-keys -t "$P1" "echo 'Starting R Console for chpt1 Analysis...'" C-m
tmux send-keys -t "$P1" "echo '==========================================='" C-m
tmux send-keys -t "$P1" "sleep 1" C-m

# Start R interactively
tmux send-keys -t "$P1" "R" C-m

# ---- Pane 2: GPU monitor (auto or chosen) ----
case "$GPU_WATCH" in
  auto)
    if have gpustat; then
      tmux send-keys -t "$P2" "watch -n1 gpustat --color" C-m
    elif have nvidia-smi; then
      tmux send-keys -t "$P2" "watch -n1 nvidia-smi" C-m
    elif have nvtop; then
      tmux send-keys -t "$P2" "nvtop" C-m
    else
      tmux send-keys -t "$P2" "echo 'No GPU tool found; showing CPU/IO overview instead'; vmstat 1" C-m
    fi
    ;;
  none)
    tmux send-keys -t "$P2" "bash" C-m
    ;;
  *)
    # Allow any custom command
    tmux send-keys -t "$P2" "$GPU_WATCH" C-m
    ;;
esac

# ---- Pane 3: working shell (with conda support) ----
cmd_p3="cd \"$WORKDIR\"; "
if [ -n "$CONDA_ENV" ]; then
  # Initialize conda in this shell if available, then activate
  cmd_p3+="if command -v conda >/dev/null 2>&1; then __conda_setup=\"\$(conda shell.bash hook 2>/dev/null)\" && eval \"\$__conda_setup\"; fi; "
  cmd_p3+="conda activate \"$CONDA_ENV\" || echo 'Could not activate: $CONDA_ENV'; "
fi
# Load modules in working shell too for convenience
cmd_p3+="module purge; "
cmd_p3+="module load rh9/R/4.5.0; "
cmd_p3+="module load rh9/gdal/3.11.0; "
cmd_p3+="module load rh9/geos 2>/dev/null || true; "
cmd_p3+="module load rh9/proj 2>/dev/null || true; "
cmd_p3+="module load rh9/udunits 2>/dev/null || true; "
cmd_p3+="export PROJECT_ROOT=/gpfs/data1/vclgp/lmaden/chpt1; "
cmd_p3+="clear; "
cmd_p3+="echo 'Working shell - modules loaded, PROJECT_ROOT set'"

tmux send-keys -t "$P3" "$cmd_p3" C-m

# Focus pane 1 (R console) and attach/switch
tmux select-pane -t "$P1"
attach_or_switch