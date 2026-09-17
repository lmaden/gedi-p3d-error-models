#!/usr/bin/env bash
# =====================================================================
# run_loo_queue.sh  --  chew through LOO folds, one at a time
#
# WHAT IT DOES
#   Works down a list of folds. Before starting one it "claims" it by
#   creating a directory, which either succeeds or fails atomically. So
#   you can run this on as many nodes as you like at the same time and no
#   two will ever pick the same fold. No scheduler, no coordination.
#
# HOW TO USE IT
#   On each node you have access to:
#
#     export PROJECT_ROOT=/gpfs/data1/vclgp/lmaden/chpt1
#     cd $PROJECT_ROOT/scripts/reviewed
#     nohup bash run_loo_queue.sh > queue_$(hostname -s).out 2>&1 &
#
#   Each instance uses 64 cores (4 chains x 16 threads). To run two folds
#   at once on one node, just launch it twice; the claim system handles it.
#   Only do that if the admin confirms 128 cores on one node is acceptable.
#
# WATCHING IT
#   bash run_loo_queue.sh --status        prints progress and exits
#   tail -n 5 $PROJECT_ROOT/models/loo/runlogs/dtm_site18_*.out
#
# IF A FOLD FAILS
#   Its claim is renamed to NAME.failed and the runner moves on rather
#   than retrying blindly. To retry after fixing whatever broke:
#     rm -rf $PROJECT_ROOT/models/loo/queue/dtm_site7.failed
#
# STOPPING
#   pkill -f run_loo_queue          stops the runner after the current
#                                    fold; it does not kill a running fit
#   pkill -f loo_fold_dtm           kills the fit itself
# =====================================================================

set -u

PROJECT_ROOT=${PROJECT_ROOT:-/gpfs/data1/vclgp/lmaden/chpt1}
export PROJECT_ROOT
HOST=$(hostname -s)

# HOME is node-local, so ~/.Renviron exists only on gsapp22. Everything it
# points at lives on GPFS and is visible everywhere, so set it explicitly.
export R_LIBS=${R_LIBS:-/gpfs/data1/vclgp/lmaden/Rlib}
export CMDSTAN=${CMDSTAN:-/gpfs/data1/vclgp/lmaden/cmdstan/cmdstan-2.37.0}

# One scratch dir per node, or simultaneous Stan compilations collide.
export TMPDIR=$PROJECT_ROOT/tmp/$HOST

# 'module' is a shell function and may not survive into a script.
if ! command -v module >/dev/null 2>&1; then
  [ -f /apps/lmod/lmod/init/bash ]   && . /apps/lmod/lmod/init/bash
  [ -f /etc/profile.d/modules.sh ]   && . /etc/profile.d/modules.sh
fi
module purge 2>/dev/null || true
module load rh9/R/4.5.0 2>/dev/null || true

RSCRIPT=$(command -v Rscript 2>/dev/null)
[ -x "$RSCRIPT" ] || RSCRIPT=/apps/rh9/R/4.5.0/bin/Rscript

SCRIPTS=$PROJECT_ROOT/scripts/reviewed
QUEUE=$PROJECT_ROOT/models/loo/queue
RUNLOGS=$PROJECT_ROOT/models/loo/runlogs

mkdir -p "$QUEUE" "$RUNLOGS" "$TMPDIR"

# Ordered longest fold first. Holding out a SMALL site leaves nearly all
# the data, so those folds are the expensive ones and should start first,
# otherwise they end up as a long tail at the end.
#   site:  18   17   14   15    2    20     3     8     13    12
#   rows:  23   27   413  488   514  1333   2314  2419  2546  3406
#   site:  9    11   4    7     6    19     5      1
#   rows:  3824 4137 4440 5988  6265 16533  71111  76631
DTM_SITES="18 17 14 15 2 20 3 8 13 12 9 11 4 7 6 19 5 1"

# ---- status mode ----------------------------------------------------
if [ "${1:-}" = "--status" ]; then
  echo "LOO queue status at $(date)"
  echo "queue dir: $QUEUE"
  printf "%-14s %s\n" "FOLD" "STATE"
  for S in $DTM_SITES; do
    T="dtm_site${S}"
    if   [ -e "$QUEUE/${T}.done"   ]; then ST="done"
    elif [ -d "$QUEUE/${T}.failed" ]; then ST="FAILED"
    elif [ -d "$QUEUE/${T}.claim"  ]; then ST="running on $(cat "$QUEUE/${T}.claim/owner" 2>/dev/null)"
    else ST="waiting"; fi
    printf "%-14s %s\n" "$T" "$ST"
  done
  echo
  echo "done:    $(ls -d "$QUEUE"/*.done   2>/dev/null | wc -l)"
  echo "running: $(ls -d "$QUEUE"/*.claim  2>/dev/null | wc -l)"
  echo "failed:  $(ls -d "$QUEUE"/*.failed 2>/dev/null | wc -l)"
  exit 0
fi

cd "$SCRIPTS" || { echo "FATAL: cannot cd to $SCRIPTS"; exit 1; }
if [ ! -f loo_fold_dtm.R ]; then
  echo "FATAL: loo_fold_dtm.R not found in $SCRIPTS"; exit 1
fi
if [ ! -x "$RSCRIPT" ]; then
  echo "FATAL: no usable Rscript. Tried PATH and /apps/rh9/R/4.5.0/bin/Rscript"; exit 1
fi

# Other people share these machines. A fold wants 64 cores; do not start one
# unless there is room, unless explicitly forced.
CORES=$(nproc --all)
LOAD=$(awk '{printf "%d", $1}' /proc/loadavg)
FREE=$(( CORES - LOAD ))
echo "capacity: $CORES cores, load $LOAD, about $FREE free"
if [ "$FREE" -lt 64 ] && [ "${1:-}" != "--force" ]; then
  echo ""
  echo "REFUSING TO START: fewer than 64 cores free on $HOST."
  echo "Other users are on this machine and a fold would oversubscribe it."
  echo "Wait for load to fall, or override with:  bash run_loo_queue.sh --force"
  exit 1
fi

echo "======================================================================"
echo "LOO queue runner"
echo "host:         $HOST"
echo "started:      $(date)"
echo "PROJECT_ROOT: $PROJECT_ROOT"
echo "TMPDIR:       $TMPDIR"
echo "R_LIBS:       $R_LIBS"
echo "CMDSTAN:      $CMDSTAN"
echo "Rscript:      $RSCRIPT"
echo "each fold:    4 chains x 16 threads = 64 cores, 3000 iterations"
echo "======================================================================"

for SITE in $DTM_SITES; do
  TASK="dtm_site${SITE}"
  DONE="$QUEUE/${TASK}.done"
  CLAIM="$QUEUE/${TASK}.claim"
  FAILED="$QUEUE/${TASK}.failed"

  [ -e "$DONE"   ] && { echo "[skip] $TASK already done";   continue; }
  [ -d "$FAILED" ] && { echo "[skip] $TASK previously failed, clear it to retry"; continue; }

  # mkdir either creates the directory or fails. That is the whole lock.
  if ! mkdir "$CLAIM" 2>/dev/null; then
    echo "[skip] $TASK claimed by another runner"
    continue
  fi
  echo "$HOST pid $$ $(date)" > "$CLAIM/owner"

  echo
  echo "---------------------------------------------------------------"
  echo "[start] $TASK on $HOST at $(date)"
  echo "---------------------------------------------------------------"
  T0=$(date +%s)

  nice -n 10 "$RSCRIPT" loo_fold_dtm.R --site "$SITE" \
      > "$RUNLOGS/${TASK}_${HOST}.out" 2>&1
  RC=$?

  MINS=$(( ( $(date +%s) - T0 ) / 60 ))
  if [ $RC -eq 0 ]; then
    touch "$DONE"
    rm -rf "$CLAIM"
    echo "[done]  $TASK in ${MINS} min"
    grep -E "^  (divergences|fold diagnostics|RMSE,|cov90,)" \
         "$RUNLOGS/${TASK}_${HOST}.out" 2>/dev/null
  else
    mv "$CLAIM" "$FAILED" 2>/dev/null
    echo "[FAIL]  $TASK exited $RC after ${MINS} min"
    echo "        see $RUNLOGS/${TASK}_${HOST}.out"
    tail -n 15 "$RUNLOGS/${TASK}_${HOST}.out"
  fi
done

echo
echo "======================================================================"
echo "runner on $HOST finished at $(date)"
echo "nothing left that this runner can claim"
echo "======================================================================"
