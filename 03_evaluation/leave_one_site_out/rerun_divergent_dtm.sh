#!/usr/bin/env bash
# =====================================================================
# rerun_divergent_dtm.sh
#
# Reruns the six DTM folds that had divergent transitions at
# adapt_delta 0.98, this time at 0.99.
#
#   site 15  (488 rows held out)   3 divergences
#   site 14  (413)                 3
#   site 13  (2545)                1
#   site 4   (4439)                1
#   site 7   (5988)                1
#   site 8   (2416)                1
#
# The original results are NOT touched. Everything this produces carries
# an _ad099 suffix, so both versions sit side by side on disk and the
# choice of which to report stays visible.
#
# ONE FOLD PER NODE, 64 cores. Same footprint as the original run.
#
# HOW TO USE
#   On each node, one at a time:
#     cd /gpfs/data1/vclgp/lmaden/chpt1/scripts/reviewed
#     nohup bash rerun_divergent_dtm.sh > rerun_$(hostname -s).out 2>&1 &
#
#   Progress:  bash rerun_divergent_dtm.sh --status
#   Stop:      pkill -f rerun_divergent_dtm     (after current fold)
#              pkill -f loo_fold_dtm            (kills the fold too)
#
# EXPECTED TIME
#   About 173 hours of work in total, since 0.99 costs roughly 25% more
#   than the 138 hours these six took at 0.98. On two nodes that is
#   roughly four days, set by the longest fold rather than the average.
# =====================================================================

set -u

PROJECT_ROOT=${PROJECT_ROOT:-/gpfs/data1/vclgp/lmaden/chpt1}
export PROJECT_ROOT
HOST=$(hostname -s)

export R_LIBS=${R_LIBS:-/gpfs/data1/vclgp/lmaden/Rlib}
export CMDSTAN=${CMDSTAN:-/gpfs/data1/vclgp/lmaden/cmdstan/cmdstan-2.37.0}
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
QUEUE=$PROJECT_ROOT/models/loo/queue_ad099
RUNLOGS=$PROJECT_ROOT/models/loo/runlogs
mkdir -p "$QUEUE" "$RUNLOGS" "$TMPDIR"

ADAPT_DELTA=0.99
# Longest fold first, so the expensive ones are not left to the end.
SITES="15 13 14 4 8 7"

# ---- status ---------------------------------------------------------
if [ "${1:-}" = "--status" ]; then
  echo "adapt_delta $ADAPT_DELTA rerun status at $(date)"
  for S in $SITES; do
    T="dtm_site${S}_ad099"
    if   [ -e "$QUEUE/${T}.done"   ]; then ST="done"
    elif [ -d "$QUEUE/${T}.failed" ]; then ST="FAILED"
    elif [ -d "$QUEUE/${T}.claim"  ]; then ST="running on $(cat "$QUEUE/${T}.claim/owner" 2>/dev/null)"
    else ST="waiting"; fi
    printf "  %-22s %s\n" "$T" "$ST"
  done
  echo "done: $(ls -d "$QUEUE"/*.done 2>/dev/null | wc -l) of 6"
  exit 0
fi

cd "$SCRIPTS" || { echo "FATAL: cannot cd to $SCRIPTS"; exit 1; }
[ -f loo_fold_dtm.R ] || { echo "FATAL: loo_fold_dtm.R not found"; exit 1; }
[ -x "$RSCRIPT" ]     || { echo "FATAL: no usable Rscript"; exit 1; }

# The patched script must support --adapt-delta, or this would silently
# rerun at 0.98 and waste four days.
if ! grep -q -- "--adapt-delta" loo_fold_dtm.R; then
  echo "FATAL: loo_fold_dtm.R does not support --adapt-delta."
  echo "Upload the patched version before running this."
  exit 1
fi

CORES=$(nproc --all)
LOAD=$(awk '{printf "%d", $1}' /proc/loadavg)
FREE=$(( CORES - LOAD ))
echo "capacity: $CORES cores, load $LOAD, about $FREE free"
if [ "$FREE" -lt 64 ] && [ "${1:-}" != "--force" ]; then
  echo "REFUSING TO START: under 64 cores free on $HOST. Retry later, or --force."
  exit 1
fi

echo "======================================================================"
echo "DTM divergent-fold rerun at adapt_delta $ADAPT_DELTA"
echo "host:     $HOST"
echo "started:  $(date)"
echo "sites:    $SITES"
echo "outputs:  *_ad099.*  (originals untouched)"
echo "======================================================================"

for SITE in $SITES; do
  TASK="dtm_site${SITE}_ad099"
  DONE="$QUEUE/${TASK}.done"; CLAIM="$QUEUE/${TASK}.claim"; FAILED="$QUEUE/${TASK}.failed"

  [ -e "$DONE"   ] && { echo "[skip] $TASK already done"; continue; }
  [ -d "$FAILED" ] && { echo "[skip] $TASK previously failed"; continue; }
  mkdir "$CLAIM" 2>/dev/null || { echo "[skip] $TASK claimed elsewhere"; continue; }
  echo "$HOST pid $$ $(date)" > "$CLAIM/owner"

  echo ""
  echo "[start] $TASK on $HOST at $(date)"
  T0=$(date +%s)
  nice -n 10 "$RSCRIPT" loo_fold_dtm.R --site "$SITE" --adapt-delta "$ADAPT_DELTA" \
      > "$RUNLOGS/${TASK}_${HOST}.out" 2>&1
  RC=$?
  MINS=$(( ( $(date +%s) - T0 ) / 60 ))

  if [ $RC -eq 0 ]; then
    touch "$DONE"; rm -rf "$CLAIM"
    echo "[done]  $TASK in ${MINS} min"
    grep -E "^  (divergences|fold diagnostics|RMSE,|cov90,)" \
         "$RUNLOGS/${TASK}_${HOST}.out" 2>/dev/null
  else
    mv "$CLAIM" "$FAILED" 2>/dev/null
    echo "[FAIL]  $TASK exited $RC after ${MINS} min"
    tail -n 15 "$RUNLOGS/${TASK}_${HOST}.out"
  fi
done

echo ""
echo "runner on $HOST finished at $(date)"
