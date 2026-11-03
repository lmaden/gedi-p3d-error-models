#!/usr/bin/env bash
# =====================================================================
# run_rework_queue.sh  --  the 21 model runs of the error models
# Paper: runs the fits behind Table 3, Table S2, Table S5, Fig. S1 and Text S2.
#
# 3 full fits (chm, dtm, dsm; rework_fullfit.R), then 3 x K folds
# (K from manuscript_tables/rework/folds.csv; rework_fold.R). Claim-directory
# locking, so several runners (on one node or several) never pick the same
# job. Each job uses 32 cores (4 chains x 8 threads).
#
# BEFORE LAUNCHING: checkpoints/09_rework_prep.rds (rework_model_prep.R).
#
# USE (one runner per 32 free cores; on a big node start two or three)
#   export PROJECT_ROOT=<CLUSTER_ROOT>/chpt1 R_LIBS=<R_LIBRARY> CMDSTAN=<CLUSTER_ROOT>/cmdstan/cmdstan-2.37.0
#   cd $PROJECT_ROOT/scripts/reviewed
#   nohup bash run_rework_queue.sh > queue_rework_$(hostname -s)_$RANDOM.out 2>&1 &
# STATUS
#   bash run_rework_queue.sh --status
# RETRY A FAILED JOB (with sampler overrides; the fit file name records them,
# so a retry never reuses the failed fit)
#   echo "--init 0" > $PROJECT_ROOT/models/rework/queue/<job>.args      (or --adapt-delta 0.995)
#   rm -rf $PROJECT_ROOT/models/rework/queue/<job>.failed   (then relaunch a runner)
# SHRINK BY ONE RUNNER (it finishes its current job first)
#   touch $PROJECT_ROOT/models/rework/queue/drain
# STOP
#   pkill -f run_rework_queue    (the runners; their fits keep running, see below)
#   pkill -f rework_f            (kills the R session of the running fit)
#   pkill -f model_6f7e908372fe057d52b855d76c4deefd   (the Stan chains survive the R
#                                 session and must be killed too; hash of this model)
# =====================================================================

set -u
PROJECT_ROOT=${PROJECT_ROOT:-<CLUSTER_ROOT>/chpt1}
export PROJECT_ROOT
HOST=$(hostname -s)
export R_LIBS=${R_LIBS:-<R_LIBRARY>}
export CMDSTAN=${CMDSTAN:-<CLUSTER_ROOT>/cmdstan/cmdstan-2.37.0}
export TMPDIR=$PROJECT_ROOT/tmp/$HOST

if ! command -v module >/dev/null 2>&1; then
  [ -f /apps/lmod/lmod/init/bash ] && . /apps/lmod/lmod/init/bash
  [ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
fi
module purge 2>/dev/null || true
module load rh9/R/4.5.0 2>/dev/null || true
RSCRIPT=$(command -v Rscript 2>/dev/null)
[ -x "$RSCRIPT" ] || RSCRIPT=/apps/rh9/R/4.5.0/bin/Rscript

SCRIPTS=$PROJECT_ROOT/scripts/reviewed
QUEUE=$PROJECT_ROOT/models/rework/queue
RUNLOGS=$PROJECT_ROOT/models/rework/runlogs
PREP=$PROJECT_ROOT/checkpoints/09_rework_prep.rds
mkdir -p "$QUEUE" "$RUNLOGS" "$TMPDIR"

# Full fits first (they carry the explained-variance results), then the folds. Fold 1 holds
# the largest site (internal id 1), so it trains on the least data and runs fastest.
FOLDS_CSV=$PROJECT_ROOT/manuscript_tables/rework/folds.csv
[ -f "$FOLDS_CSV" ] || { echo "FATAL: $FOLDS_CSV missing. Run rework_model_prep.R first."; exit 1; }
KS=$(tail -n +2 "$FOLDS_CSV" | cut -d, -f3 | sort -n -u)
# fold 1 holds the largest sites, so it trains on the least data: run it last
KORDER="$(echo "$KS" | grep -vx 1 | tr '\n' ' ') 1"
JOBS="fit_dtm fit_dsm fit_chm"
for K in $KORDER; do
  for R in dtm dsm chm; do JOBS="$JOBS ${R}_fold${K}"; done
done
NJOBS=$(echo $JOBS | wc -w)

if [ "${1:-}" = "--status" ]; then
  echo "rework queue at $(date)"
  for J in $JOBS; do
    if   [ -e "$QUEUE/${J}.done"   ]; then ST="done"
    elif [ -d "$QUEUE/${J}.failed" ]; then ST="FAILED"
    elif [ -d "$QUEUE/${J}.claim"  ]; then ST="running on $(cat "$QUEUE/${J}.claim/owner" 2>/dev/null)"
    else ST="waiting"; fi
    printf "%-12s %s\n" "$J" "$ST"
  done
  echo "done $(ls -d "$QUEUE"/*.done 2>/dev/null | wc -l) | running $(ls -d "$QUEUE"/*.claim 2>/dev/null | wc -l) | failed $(ls -d "$QUEUE"/*.failed 2>/dev/null | wc -l) | of $NJOBS"
  exit 0
fi

cd "$SCRIPTS" || { echo "FATAL: cannot cd to $SCRIPTS"; exit 1; }
for F in rework_fullfit.R rework_fold.R; do
  [ -f "$F" ] || { echo "FATAL: $F not in $SCRIPTS"; exit 1; }
done
[ -f "$PREP" ] || { echo "FATAL: $PREP missing. Run rework_model_prep.R first."; exit 1; }

CORES=$(nproc --all); LOAD=$(awk '{printf "%d", $1}' /proc/loadavg); FREE=$(( CORES - LOAD ))
echo "capacity: $CORES cores, load $LOAD, about $FREE free"
if [ "$FREE" -lt 32 ] && [ "${1:-}" != "--force" ]; then
  echo "REFUSING TO START: fewer than 32 cores free on $HOST (override: --force)"; exit 1
fi
echo "rework queue runner on $HOST, started $(date)"

for J in $JOBS; do
  # drain: `touch $QUEUE/drain` retires ONE runner after its current job (rm is atomic)
  rm "$QUEUE/drain" 2>/dev/null && { echo "[drain] runner pid $$ on $HOST retiring at $(date)"; exit 0; }
  DONE="$QUEUE/${J}.done"; CLAIM="$QUEUE/${J}.claim"; FAILED="$QUEUE/${J}.failed"
  [ -e "$DONE" ] && continue
  [ -d "$FAILED" ] && { echo "[skip] $J failed earlier"; continue; }
  mkdir "$CLAIM" 2>/dev/null || continue
  echo "$HOST pid $$ $(date)" > "$CLAIM/owner"
  EXTRA=""; [ -f "$QUEUE/${J}.args" ] && EXTRA=$(cat "$QUEUE/${J}.args")
  echo "[start] $J on $HOST at $(date) ${EXTRA:+(overrides: $EXTRA)}"
  T0=$(date +%s)
  case $J in
    fit_*) R=${J#fit_}; nice -n 10 "$RSCRIPT" rework_fullfit.R --response "$R" $EXTRA > "$RUNLOGS/${J}_${HOST}.out" 2>&1 ;;
    *)     R=${J%_fold*}; K=${J##*_fold}
           nice -n 10 "$RSCRIPT" rework_fold.R --response "$R" --fold "$K" $EXTRA > "$RUNLOGS/${J}_${HOST}.out" 2>&1 ;;
  esac
  RC=$?
  MINS=$(( ( $(date +%s) - T0 ) / 60 ))
  if [ $RC -eq 0 ]; then
    touch "$DONE"; rm -rf "$CLAIM"; echo "[done]  $J in ${MINS} min"
    grep -E "FULL FIT DIAGNOSTICS|FOLD DIAGNOSTICS|model .* vs constant|R2 total" "$RUNLOGS/${J}_${HOST}.out" 2>/dev/null
  else
    mv "$CLAIM" "$FAILED" 2>/dev/null
    echo "[FAIL]  $J exited $RC after ${MINS} min; see $RUNLOGS/${J}_${HOST}.out"
    tail -n 12 "$RUNLOGS/${J}_${HOST}.out"
  fi
done
echo "runner on $HOST finished at $(date)"
