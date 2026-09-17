#!/usr/bin/env bash
# =====================================================================
# als_pointcloud_inventory.sh
#
# v2 of the ALS point-cloud inventory. v1 used dir-name patterns
# (-iname "laz", -iname "*point_cloud*", etc.) to find candidate
# subdirs, which missed LAZ_ground because the exact-name "laz" pattern
# does not match "LAZ_ground". v2 instead searches by file content:
# locates all *.las / *.laz / *.copc.laz files under each site root,
# groups by parent directory, and reports counts + REAL disk sizes
# (via GNU find -printf "%s") + sample filenames.
#
# This handles arbitrary directory naming and also surfaces the case
# where files are stubs / symlinks / index entries (small total size
# despite large file count).
#
# Output: one row per (site, point_cloud_dir). A site with no PC files
# gets one row with empty pc_dir. A site with multiple PC dirs gets
# multiple rows.
#
# Usage:
#   bash als_pointcloud_inventory.sh
#
# Outputs:
#   $PROJECT_ROOT/manuscript_tables/groundwork_task2_pointcloud_inventory_v2.csv
#   $PROJECT_ROOT/manuscript_tables/groundwork_task2_pointcloud_inventory_v2.log
# =====================================================================

set -uo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/gpfs/data1/vclgp/lmaden/chpt1}"
OUT_DIR="${PROJECT_ROOT}/manuscript_tables"
OUT_CSV="${OUT_DIR}/groundwork_task2_pointcloud_inventory_v2.csv"
LOG="${OUT_DIR}/groundwork_task2_pointcloud_inventory_v2.log"

mkdir -p "$OUT_DIR"

# Per-site root paths (parents of the /chm subdirs from tracker.xlsx).
declare -A SITE_ROOT
SITE_ROOT[1]="/gpfs/data1/vclgp/data/gedi/imported/usa/usda_me"
SITE_ROOT[2]="/gpfs/data1/vclgp/data/gedi/imported/usa/nasa_howland"
SITE_ROOT[3]="/gpfs/data1/vclgp/data/gedi/imported/usa/neon_sawb"
SITE_ROOT[4]="/gpfs/data1/vclgp/data/gedi/imported/usa/neon_harv2019"
SITE_ROOT[5]="/gpfs/data1/vclgp/data/gedi/imported/usa/usda_sc"
SITE_ROOT[6]="/gpfs/data1/vclgp/data/gedi/imported/usa/neon_jerc2021"
SITE_ROOT[7]="/gpfs/data1/vclgp/data/gedi/imported/usa/neon_tall2021"
SITE_ROOT[8]="/gpfs/data1/vclgp/data/gedi/imported/usa/neon_dela2021"
SITE_ROOT[9]="/gpfs/data1/vclgp/data/gedi/imported/usa/neon_leno2021"
SITE_ROOT[10]="/gpfs/data1/vclgp/data/gedi/imported/usa/neon_clbj2022"
SITE_ROOT[11]="/gpfs/data1/vclgp/data/gedi/imported/usa/neon_konz2020"
SITE_ROOT[12]="/gpfs/data1/vclgp/data/gedi/imported/usa/neon_stei2022"
SITE_ROOT[13]="/gpfs/data1/vclgp/data/gedi/imported/usa/neon_unde2022"
SITE_ROOT[14]="/gpfs/data1/vclgp/data/gedi/imported/usa/neon_steicheq2022"
SITE_ROOT[15]="/gpfs/data1/vclgp/data/gedi/imported/usa/neon_wood2021"
SITE_ROOT[17]="/gpfs/data1/vclgp/data/gedi/imported/usa/neon_ster2022"
SITE_ROOT[18]="/gpfs/data1/vclgp/data/gedi/imported/usa/neon_cper2021"
SITE_ROOT[19]="/gpfs/data1/vclgp/data/gedi/imported/usa/ltbmu_20180710"
SITE_ROOT[20]="/gpfs/data1/vclgp/data/gedi/imported/usa/neon_sjer"

# Manuscript site -> tracker site
declare -A MS_TO_TR
MS_TO_TR[1]=1;  MS_TO_TR[2]=2;  MS_TO_TR[3]=3;  MS_TO_TR[4]=4
MS_TO_TR[5]=5;  MS_TO_TR[6]=6;  MS_TO_TR[7]=7;  MS_TO_TR[8]=8
MS_TO_TR[9]=9;  MS_TO_TR[10]=10; MS_TO_TR[11]=11; MS_TO_TR[12]=12
MS_TO_TR[13]=13; MS_TO_TR[14]=14; MS_TO_TR[15]=15
MS_TO_TR[16]=17; MS_TO_TR[17]=18; MS_TO_TR[18]=19; MS_TO_TR[19]=20

# Manuscript site -> short name
declare -A MS_TO_SHORT
MS_TO_SHORT[1]="usda_me"
MS_TO_SHORT[2]="nasa_howland"
MS_TO_SHORT[3]="neon_sawb"
MS_TO_SHORT[4]="neon_harv2019"
MS_TO_SHORT[5]="usda_sc"
MS_TO_SHORT[6]="neon_jerc2021"
MS_TO_SHORT[7]="neon_tall2021"
MS_TO_SHORT[8]="neon_dela2021"
MS_TO_SHORT[9]="neon_leno2021"
MS_TO_SHORT[10]="neon_clbj2022"
MS_TO_SHORT[11]="neon_konz2020"
MS_TO_SHORT[12]="neon_stei2022"
MS_TO_SHORT[13]="neon_unde2022"
MS_TO_SHORT[14]="neon_steicheq2022"
MS_TO_SHORT[15]="neon_wood2021"
MS_TO_SHORT[16]="neon_ster2022"
MS_TO_SHORT[17]="neon_cper2021"
MS_TO_SHORT[18]="ltbmu_20180710"
MS_TO_SHORT[19]="neon_sjer"

# Site flag status
declare -A MS_FLAG
MS_FLAG[1]="FLAGGED"
MS_FLAG[2]="FLAGGED"
MS_FLAG[3]="FLAGGED"
for ms in 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; do MS_FLAG[$ms]="OK"; done

ts() { date '+%Y-%m-%d %H:%M:%S'; }

: > "$LOG"
echo "[$(ts)] === ALS point-cloud inventory (file-content based) ===" | tee -a "$LOG"
echo "[$(ts)] PROJECT_ROOT = $PROJECT_ROOT" | tee -a "$LOG"
echo "[$(ts)] OUT_CSV      = $OUT_CSV" | tee -a "$LOG"
echo "" | tee -a "$LOG"

# CSV header
echo "manuscript_site,tracker_site,site_short_name,flag_status,site_root,site_root_exists,pc_dir,n_files,n_las,n_laz,n_copc,total_gb,n_symlinks,sample_1,sample_2,sample_3" > "$OUT_CSV"

inventory_one_site() {
  local ms=$1
  local tr="${MS_TO_TR[$ms]}"
  local short="${MS_TO_SHORT[$ms]}"
  local flag="${MS_FLAG[$ms]}"
  local root="${SITE_ROOT[$tr]}"

  if [[ ! -d "$root" ]]; then
    echo "  [ms $ms / tr $tr] $short ($flag): SITE ROOT NOT FOUND ($root)" | tee -a "$LOG"
    echo "$ms,$tr,$short,$flag,$root,0,,0,0,0,0,0,0,,," >> "$OUT_CSV"
    return
  fi

  # Find all .las/.laz/.copc.laz files under root, maxdepth 4.
  # File-content search; directory naming convention does not matter.
  local all_files
  all_files=$(find "$root" -maxdepth 4 -type f \
    \( -iname "*.las" -o -iname "*.laz" -o -iname "*.copc.laz" \) \
    2>/dev/null)

  if [[ -z "$all_files" ]]; then
    # Fallback: also check for symlinks
    all_files=$(find "$root" -maxdepth 4 -type l \
      \( -iname "*.las" -o -iname "*.laz" -o -iname "*.copc.laz" \) \
      2>/dev/null)
  fi

  if [[ -z "$all_files" ]]; then
    echo "  [ms $ms / tr $tr] $short ($flag): NO point-cloud files found under $root" | tee -a "$LOG"
    echo "$ms,$tr,$short,$flag,$root,1,,0,0,0,0,0,0,,," >> "$OUT_CSV"
    return
  fi

  # Group by parent directory
  local pc_dirs
  pc_dirs=$(echo "$all_files" | xargs -L 1 dirname 2>/dev/null | sort -u)

  while IFS= read -r d; do
    [[ -z "$d" ]] && continue

    local n_files n_las n_laz n_copc n_symlinks total_bytes total_gb
    local sample_1="" sample_2="" sample_3=""

    # Counts (regular files + symlinks both counted, since GPFS
    # commonly uses symlinks for staged data).
    n_files=$(find "$d" -maxdepth 1 \( -type f -o -type l \) \
      \( -iname "*.las" -o -iname "*.laz" -o -iname "*.copc.laz" \) \
      2>/dev/null | wc -l)
    n_las=$(find "$d" -maxdepth 1 \( -type f -o -type l \) -iname "*.las" \
      2>/dev/null | wc -l)
    n_laz=$(find "$d" -maxdepth 1 \( -type f -o -type l \) -iname "*.laz" \
      2>/dev/null | wc -l)
    n_copc=$(find "$d" -maxdepth 1 \( -type f -o -type l \) -iname "*.copc.laz" \
      2>/dev/null | wc -l)
    n_symlinks=$(find "$d" -maxdepth 1 -type l \
      \( -iname "*.las" -o -iname "*.laz" -o -iname "*.copc.laz" \) \
      2>/dev/null | wc -l)

    # Real disk size, following symlinks via -L so we get target size
    # not link size. -printf "%s\n" reports byte count of the target.
    total_bytes=$(find -L "$d" -maxdepth 1 -type f \
      \( -iname "*.las" -o -iname "*.laz" -o -iname "*.copc.laz" \) \
      -printf "%s\n" 2>/dev/null | awk '{ sum += $1 } END { print sum+0 }')
    total_gb=$(awk -v b="$total_bytes" 'BEGIN { printf "%.4f", b / (1024^3) }')

    # Sample 3 filenames
    local samples
    samples=$(find "$d" -maxdepth 1 \( -type f -o -type l \) \
      \( -iname "*.las" -o -iname "*.laz" -o -iname "*.copc.laz" \) \
      2>/dev/null | sort | head -3)
    sample_1=$(echo "$samples" | sed -n '1p' | xargs -L 1 basename 2>/dev/null)
    sample_2=$(echo "$samples" | sed -n '2p' | xargs -L 1 basename 2>/dev/null)
    sample_3=$(echo "$samples" | sed -n '3p' | xargs -L 1 basename 2>/dev/null)

    printf "  [ms %2d / tr %2d] %-20s (%s): pc_dir=%s\n" "$ms" "$tr" "$short" "$flag" "$d" | tee -a "$LOG"
    printf "      n=%d  las=%d  laz=%d  copc=%d  symlinks=%d  size=%s GB\n" \
      "$n_files" "$n_las" "$n_laz" "$n_copc" "$n_symlinks" "$total_gb" | tee -a "$LOG"
    printf "      sample: %s\n" "$sample_1" | tee -a "$LOG"

    # CSV row (quote pc_dir and samples to handle any spaces)
    echo "$ms,$tr,$short,$flag,$root,1,\"$d\",$n_files,$n_las,$n_laz,$n_copc,$total_gb,$n_symlinks,\"$sample_1\",\"$sample_2\",\"$sample_3\"" >> "$OUT_CSV"
  done <<< "$pc_dirs"
}

for ms in $(seq 1 19); do
  inventory_one_site "$ms"
done

echo "" | tee -a "$LOG"
echo "[$(ts)] === Inventory complete ===" | tee -a "$LOG"
echo "" | tee -a "$LOG"

echo "Summary by site (manuscript number, total files across all PC dirs, total GB):" | tee -a "$LOG"
awk -F, 'NR>1 {
  ms = $1; short = $3; flag = $4
  n = $8 + 0; gb = $12 + 0
  total_n[ms"\t"short"\t"flag] += n
  total_gb[ms"\t"short"\t"flag] += gb
  ms_seen[ms"\t"short"\t"flag] = 1
}
END {
  for (k in ms_seen) {
    printf "  %s\t  n=%d  total_gb=%.3f\n", k, total_n[k], total_gb[k]
  }
}' "$OUT_CSV" | sort -k1,1n | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "Diagnostic checks to run after this:" | tee -a "$LOG"
echo "  1. If total_gb is suspiciously small (<1 GB per site), check whether files are stubs:" | tee -a "$LOG"
echo "     ls -lL <one_file.laz>  # follow symlinks; shows real target size" | tee -a "$LOG"
echo "     file <one_file.laz>    # file type (LAS/LAZ binary vs ASCII text)" | tee -a "$LOG"
echo "  2. If LAStools available, verify a file is a real point cloud:" | tee -a "$LOG"
echo "     lasinfo <one_file.laz> | head -30" | tee -a "$LOG"
echo "  3. Check classification distribution to confirm full point cloud (not just ground):" | tee -a "$LOG"
echo "     lasinfo -i <one_file.laz> -no_check 2>&1 | grep -A 20 'classification'" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "[$(ts)] === Done ===" | tee -a "$LOG"
