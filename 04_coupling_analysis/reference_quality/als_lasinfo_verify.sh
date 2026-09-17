#!/usr/bin/env bash
# =====================================================================
# als_lasinfo_verify.sh
#
# v2: detects the LAStools binary name automatically. v1 assumed
# `lasinfo` would be on PATH after `module load rh9/lastools/250710`,
# but recent LAStools packaging uses `lasinfo64` (64-bit suffix). v2
# searches for `lasinfo`, `lasinfo64`, and `lasinfo32` in order and
# uses whichever is found.
#
# Verify ALS point-cloud format on one sample file per (site, pc_dir).
# Reads groundwork_task2_pointcloud_inventory_v2.csv, runs lasinfo on
# the first sample file in each row, parses out:
#   - total point count
#   - return number distribution
#   - classification histogram
#   - bounding box
# Writes a per-row verdict CSV plus full lasinfo output log.
#
# Verdict logic:
#   - FULL_POINT_CLOUD : >= 2 distinct return numbers AND >= 2 classes
#   - LIKELY_FULL      : >= 2 returns, only 1 class (e.g. all class 1)
#   - GROUND_ONLY      : only class 2 present, 1 return
#   - LIMITED          : less than full, more than ground-only
#   - LASINFO_FAILED   : lasinfo command returned an error
#   - NOT_FOUND        : sample file path does not exist on disk
#
# Usage:
#   bash groundwork_task2_lasinfo_verify.sh
#
# Requires: lasinfo on PATH (LAStools). Per cluster setup, may need:
#   module load lastools
#
# Outputs:
#   $PROJECT_ROOT/manuscript_tables/groundwork_task2_lasinfo_verify.csv
#   $PROJECT_ROOT/manuscript_tables/groundwork_task2_lasinfo_verify.log
# =====================================================================

set -uo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/gpfs/data1/vclgp/lmaden/chpt1}"
OUT_DIR="$PROJECT_ROOT/manuscript_tables"
INPUT_CSV="$OUT_DIR/groundwork_task2_pointcloud_inventory_v2.csv"
OUT_CSV="$OUT_DIR/groundwork_task2_lasinfo_verify.csv"
LOG="$OUT_DIR/groundwork_task2_lasinfo_verify.log"

if [[ ! -f "$INPUT_CSV" ]]; then
  echo "ERROR: v2 inventory CSV not found: $INPUT_CSV"
  echo "Run als_pointcloud_inventory.sh first."
  exit 1
fi

if ! command -v lasinfo >/dev/null 2>&1 && \
   ! command -v lasinfo64 >/dev/null 2>&1 && \
   ! command -v lasinfo32 >/dev/null 2>&1; then
  echo "ERROR: no lasinfo binary found in PATH (tried lasinfo, lasinfo64, lasinfo32)"
  echo "On cluster: module load rh9/lastools/250710"
  exit 1
fi

# Pick the first available binary (LAStools 64-bit packaging uses
# lasinfo64; older builds use lasinfo).
LASINFO_BIN=""
for cand in lasinfo lasinfo64 lasinfo32; do
  if command -v "$cand" >/dev/null 2>&1; then
    LASINFO_BIN=$(command -v "$cand")
    break
  fi
done
echo "Using LASInfo binary: $LASINFO_BIN"

ts() { date '+%Y-%m-%d %H:%M:%S'; }

: > "$LOG"
echo "[$(ts)] === lasinfo verification ===" | tee -a "$LOG"
echo "[$(ts)] INPUT_CSV = $INPUT_CSV" | tee -a "$LOG"
echo "[$(ts)] OUT_CSV   = $OUT_CSV" | tee -a "$LOG"
echo "" | tee -a "$LOG"

# CSV header
echo "manuscript_site,site_short_name,flag_status,pc_dir,sample_file,n_points,n_returns,returns_present,n_classes,classes_present,verdict" > "$OUT_CSV"

verify_one() {
  local ms="$1" short="$2" flag="$3" pcdir="$4" sample="$5"

  local filepath="$pcdir/$sample"

  echo "=========================================" | tee -a "$LOG"
  echo "[ms $ms / $short ($flag)] $sample" | tee -a "$LOG"
  echo "  pc_dir: $pcdir" | tee -a "$LOG"

  if [[ ! -f "$filepath" ]]; then
    echo "  STATUS: file not found at expected path" | tee -a "$LOG"
    echo "$ms,$short,$flag,\"$pcdir\",\"$sample\",0,0,,0,,NOT_FOUND" >> "$OUT_CSV"
    return
  fi

  # Run lasinfo, capture into a temp file
  local tmpfile
  tmpfile=$(mktemp)
  if ! "$LASINFO_BIN" "$filepath" -no_check > "$tmpfile" 2>&1; then
    echo "  STATUS: lasinfo command failed" | tee -a "$LOG"
    head -30 "$tmpfile" | sed 's/^/    /' | tee -a "$LOG"
    echo "$ms,$short,$flag,\"$pcdir\",\"$sample\",0,0,,0,,LASINFO_FAILED" >> "$OUT_CSV"
    rm -f "$tmpfile"
    return
  fi

  # Parse total point count
  local n_points
  n_points=$(grep -m1 -E "number of point records" "$tmpfile" \
    | grep -oE '[0-9]+' | head -1)
  n_points=${n_points:-0}

  # Parse return number distribution
  # Format: "  number of points by return: 5678901 4321098 1234567 234567 12345"
  local returns_line returns_counts returns_present n_returns
  returns_line=$(grep -m1 "number of points by return" "$tmpfile" || true)
  returns_counts=$(echo "$returns_line" | sed -E 's/.*return:\s*//')
  returns_present=""
  n_returns=0
  if [[ -n "$returns_counts" ]]; then
    local i=1
    for cnt in $returns_counts; do
      if [[ "$cnt" -gt 0 ]] 2>/dev/null; then
        returns_present="${returns_present},${i}"
        n_returns=$((n_returns + 1))
      fi
      i=$((i + 1))
    done
    returns_present="${returns_present#,}"
  fi

  # Parse classification histogram
  # Format: lines after "histogram of classification" like
  #   "  5678901 unclassified (1)"
  #   "  4321098 ground (2)"
  local classes_present n_classes class_lines
  class_lines=$(awk '/histogram of classification/{flag=1; next}
                     flag && /^[[:space:]]*[0-9]+/{print; next}
                     flag && /^[^[:space:]]/{flag=0}' "$tmpfile" 2>/dev/null)
  classes_present=$(echo "$class_lines" | grep -oE '\([0-9]+\)' \
    | tr -d '()' | sort -nu | tr '\n' ',' | sed 's/,$//')
  if [[ -z "$classes_present" ]]; then
    # Some lasinfo versions print "X point records of class Y"
    classes_present=$(grep -oE 'class [0-9]+' "$tmpfile" \
      | grep -oE '[0-9]+' | sort -nu | tr '\n' ',' | sed 's/,$//')
  fi
  n_classes=0
  if [[ -n "$classes_present" ]]; then
    n_classes=$(echo "$classes_present" | tr ',' '\n' | grep -cv '^$')
  fi

  # Verdict
  local verdict="LIMITED"
  if [[ "$n_returns" -ge 2 && "$n_classes" -ge 2 ]]; then
    verdict="FULL_POINT_CLOUD"
  elif [[ "$n_returns" -ge 2 && "$n_classes" -le 1 ]]; then
    verdict="LIKELY_FULL"
  elif [[ "$classes_present" == "2" && "$n_returns" -le 1 ]]; then
    verdict="GROUND_ONLY"
  fi

  echo "  n_points: $n_points" | tee -a "$LOG"
  echo "  returns_present: [$returns_present]  (n=$n_returns)" | tee -a "$LOG"
  echo "  classes_present: [$classes_present]  (n=$n_classes)" | tee -a "$LOG"
  echo "  VERDICT: $verdict" | tee -a "$LOG"
  if [[ -n "$class_lines" ]]; then
    echo "  classification histogram:" | tee -a "$LOG"
    echo "$class_lines" | sed 's/^/    /' | tee -a "$LOG"
  fi

  echo "$ms,$short,$flag,\"$pcdir\",\"$sample\",$n_points,$n_returns,\"$returns_present\",$n_classes,\"$classes_present\",$verdict" >> "$OUT_CSV"

  rm -f "$tmpfile"
}

# Iterate v2 inventory CSV. Skip header and rows with no pc_dir.
# Column order in v2 CSV:
#   1=manuscript_site 2=tracker_site 3=site_short_name 4=flag_status
#   5=site_root 6=site_root_exists 7=pc_dir 8=n_files 9=n_las
#   10=n_laz 11=n_copc 12=total_gb 13=n_symlinks
#   14=sample_1 15=sample_2 16=sample_3
while IFS=',' read -r ms tr short flag root exists pcdir n nlas nlaz ncopc gb nlinks s1 s2 s3; do
  [[ "$ms" == "manuscript_site" ]] && continue

  # Strip enclosing quotes from CSV fields
  pcdir=$(echo "$pcdir" | sed 's/^"//; s/"$//')
  s1=$(echo "$s1" | sed 's/^"//; s/"$//')
  short=$(echo "$short" | sed 's/^"//; s/"$//')
  flag=$(echo "$flag" | sed 's/^"//; s/"$//')

  [[ -z "$pcdir" || -z "$s1" ]] && continue

  verify_one "$ms" "$short" "$flag" "$pcdir" "$s1"
done < "$INPUT_CSV"

echo "" | tee -a "$LOG"
echo "[$(ts)] === Verification complete ===" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "Verdict summary:" | tee -a "$LOG"
awk -F',' 'NR>1 { v=$11; gsub(/"/, "", v); count[v]++ }
           END { for (v in count) printf "  %-20s %d\n", v, count[v] }' \
  "$OUT_CSV" | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "Per-site verdicts (full table is in CSV):" | tee -a "$LOG"
awk -F',' 'NR>1 {
  ms=$1; short=$2; flag=$3; v=$11
  gsub(/"/, "", short); gsub(/"/, "", flag); gsub(/"/, "", v)
  printf "  ms %2s  %-20s  %-10s  %s\n", ms, short, flag, v
}' "$OUT_CSV" | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "[$(ts)] === Done ===" | tee -a "$LOG"
