#!/usr/bin/env bash
# =============================================================================
# als_repeat_acquisition_inventory.sh
#
# Repeat-acquisition discovery.
# Enumerate multi-year ALS CHM coverage across the 19 manuscript study sites
# on the GEDI cluster. This is the prerequisite to designing the actual
# ALS-vs-ALS comparison.
#
# Run from Pane 3 of tmux-quad-chpt1-v2.sh (gdalinfo from the loaded module
# stack is required; LAStools is not used at this stage).
#
# Outputs (under $PROJECT_ROOT/manuscript_tables/):
#   - als_repeat_acquisition_inventory.csv  -- one row per (site, sibling dir)
#   - als_repeat_acquisition_inventory.log  -- progress log
#
# Conventions:
#   - "manuscript_site" is the renumbered site (1..19) used in the manuscript.
#   - "site_dir_used_in_pipeline" is the directory name the Chapter 1 pipeline
#     points at (the canonical/known acquisition).
#   - "sibling_dir" is any directory under $ALS_BASE whose name matches the
#     site stem; rows where sibling_dir != pipeline dir are candidate
#     additional acquisitions for the multi-year comparison.
#
# Date hints are heuristic (regex over path/filename). Validate against
# provider metadata before trusting for any growth-correction calculation.
# =============================================================================

set -u  # catch unset variables, but do NOT exit on errors (no -e, no pipefail)

# Enable bash trace for debugging: run with DEBUG_TRACE=1 to see every command
if [ "${DEBUG_TRACE:-0}" = "1" ]; then
  set -x
fi

# Trap: on any error, log the line number (informational only — script continues)
trap 'echo "[TRAP] Error at line $LINENO (exit $?), continuing..." >> "${OUT_LOG:-/dev/stderr}"' ERR

PROJECT_ROOT="${PROJECT_ROOT:-/gpfs/data1/vclgp/lmaden/chpt1}"
ALS_BASE="${ALS_BASE:-/gpfs/data1/vclgp/data/gedi/imported/usa}"
OUT_DIR="${PROJECT_ROOT}/manuscript_tables"
mkdir -p "${OUT_DIR}"
OUT_CSV="${OUT_DIR}/als_repeat_acquisition_inventory.csv"
OUT_LOG="${OUT_DIR}/als_repeat_acquisition_inventory.log"

# Manuscript sites: <manuscript_num>|<site_dirname_in_pipeline>|<sibling_glob_stem>
# The stem is what we glob on; for NEON sites it's the dir name with the
# trailing year stripped, so sibling acquisitions in other years are caught.
SITES=(
  "1|usda_me|usda_me"
  "2|nasa_howland|nasa_howland"
  "3|neon_sawb|neon_sawb"
  "4|neon_harv2019|neon_harv"
  "5|usda_sc|usda_sc"
  "6|neon_jerc2021|neon_jerc"
  "7|neon_tall2021|neon_tall"
  "8|neon_dela2021|neon_dela"
  "9|neon_leno2021|neon_leno"
  "10|neon_clbj2022|neon_clbj"
  "11|neon_konz2020|neon_konz"
  "12|neon_stei2022|neon_stei"
  "13|neon_unde2022|neon_unde"
  "14|neon_steicheq2022|neon_steicheq"
  "15|neon_wood2021|neon_wood"
  "16|neon_ster2022|neon_ster"
  "17|neon_cper2021|neon_cper"
  "18|ltbmu_20180710|ltbmu"
  "19|neon_sjer|neon_sjer"
)

# CSV header
echo "manuscript_site,site_dir_used_in_pipeline,sibling_dir,is_pipeline_dir,chm_subdir_present,n_tif,total_size_mb,representative_tif,res_x_m,res_y_m,crs_epsg,extent_xmin,extent_ymin,extent_xmax,extent_ymax,date_hint" > "${OUT_CSV}"

: > "${OUT_LOG}"
log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "${OUT_LOG}"; }

# Heuristic date extraction from a string (path or filename).
extract_date_hint() {
  local s="$1"
  local d
  # YYYY-MM-DD
  d=$(echo "$s" | grep -oE '20[0-2][0-9]-[01][0-9]-[0-3][0-9]' | head -n1 || true)
  if [ -n "$d" ]; then echo "$d"; return; fi
  # YYYYMMDD (8 digits starting with 20)
  d=$(echo "$s" | grep -oE '20[0-2][0-9][01][0-9][0-3][0-9]' | head -n1 || true)
  if [ -n "$d" ]; then echo "$d"; return; fi
  # YYYY only (least specific — last resort)
  d=$(echo "$s" | grep -oE '20[0-2][0-9]' | head -n1 || true)
  echo "${d:-}"
}

probe_dir() {
  local manuscript_num="$1"
  local pipeline_dirname="$2"
  local sibling_dir="$3"
  local chm_dir="${ALS_BASE}/${sibling_dir}/chm"
  local is_pipeline; is_pipeline=$([ "$sibling_dir" = "$pipeline_dirname" ] && echo "yes" || echo "no")
  local chm_present="no"
  local n_tif=0
  local size_mb=0
  local rep_tif=""
  local res_x="" res_y="" crs="" xmin="" ymin="" xmax="" ymax=""
  local date_hint
  date_hint=$(extract_date_hint "$sibling_dir")

  if [ -d "$chm_dir" ]; then
    chm_present="yes"
    n_tif=$(find "$chm_dir" -maxdepth 1 -type f \( -iname "*.tif" -o -iname "*.tiff" \) 2>/dev/null | wc -l | tr -d ' ')
    if [ "$n_tif" -gt 0 ]; then
      size_mb=$(du -sm "$chm_dir" 2>/dev/null | awk '{print $1}' || echo "0")
      rep_tif=$(find "$chm_dir" -maxdepth 1 -type f \( -iname "*.tif" -o -iname "*.tiff" \) 2>/dev/null | head -n1)
      if [ -n "$rep_tif" ] && command -v gdalinfo >/dev/null 2>&1; then
        local info
        info=$(gdalinfo -nogcp -nomd -norat -noct "$rep_tif" 2>/dev/null || true)
        # Pixel size: "Pixel Size = (1.000000000000000,-1.000000000000000)"
        res_x=$(echo "$info" | sed -n 's/.*Pixel Size = (\([^,]*\),.*/\1/p' | head -n1)
        res_y=$(echo "$info" | sed -n 's/.*Pixel Size = ([^,]*,\([^)]*\)).*/\1/p' | head -n1)
        # EPSG: try WKT2 format first (GDAL ≥3.x default): ID["EPSG",NNNN]
        # Fall back to WKT1 format: AUTHORITY["EPSG","NNNN"]
        crs=$(echo "$info" | grep -oE 'ID\["EPSG",[0-9]+\]' | tail -n1 | grep -oE '[0-9]+' || true)
        if [ -z "$crs" ]; then
          crs=$(echo "$info" | grep -oE 'AUTHORITY\["EPSG","[0-9]+"\]' | tail -n1 | grep -oE '[0-9]+' || true)
        fi
        # Corner coords: gdalinfo prints e.g. "Lower Left  (  394420.000, 5063420.000) (...)"
        # After stripping (), comma, and converting to whitespace fields, x is $3, y is $4.
        xmin=$(echo "$info" | awk '/Lower Left/  {gsub(/[(),]/," "); print $3; exit}')
        ymin=$(echo "$info" | awk '/Lower Left/  {gsub(/[(),]/," "); print $4; exit}')
        xmax=$(echo "$info" | awk '/Upper Right/ {gsub(/[(),]/," "); print $3; exit}')
        ymax=$(echo "$info" | awk '/Upper Right/ {gsub(/[(),]/," "); print $4; exit}')
        # Prefer date hint from the tif filename if more specific than dir-name
        local fhint; fhint=$(extract_date_hint "$(basename "$rep_tif")")
        if [ -n "$fhint" ]; then date_hint="$fhint"; fi
      fi
    fi
  fi

  local rep_tif_basename=""
  [ -n "$rep_tif" ] && rep_tif_basename=$(basename "$rep_tif")

  echo "${manuscript_num},${pipeline_dirname},${sibling_dir},${is_pipeline},${chm_present},${n_tif},${size_mb},${rep_tif_basename},${res_x},${res_y},${crs},${xmin},${ymin},${xmax},${ymax},${date_hint}" >> "${OUT_CSV}"
  log "  site ${manuscript_num} (${pipeline_dirname}) <-> sibling ${sibling_dir}: chm=${chm_present} n_tif=${n_tif} crs=${crs} date_hint=${date_hint}"
}

# =============================================================================
# Main
# =============================================================================
log "=== ALS multi-year inventory ==="
log "PROJECT_ROOT = ${PROJECT_ROOT}"
log "ALS_BASE     = ${ALS_BASE}"
log "OUT_CSV      = ${OUT_CSV}"
log ""

if ! command -v gdalinfo >/dev/null 2>&1; then
  log "WARNING: gdalinfo not on PATH. Resolution/CRS/extent fields will be blank."
  log "         Source the tmux-quad-chpt1-v2.sh module stack and rerun."
fi

for entry in "${SITES[@]}"; do
  IFS='|' read -r manuscript_num pipeline_dirname stem <<< "$entry"
  log "Site ${manuscript_num} (pipeline dir: ${pipeline_dirname}, glob: ${stem}*)"

  # Collect sibling directories (while-read loop avoids mapfile portability issues)
  siblings=()
  while IFS= read -r dirpath; do
    [ -n "$dirpath" ] && siblings+=("$dirpath")
  done < <(find "${ALS_BASE}" -maxdepth 1 -type d -name "${stem}*" 2>/dev/null | sort)

  if [ "${#siblings[@]}" -eq 0 ]; then
    log "  (no directories matched ${stem}*)"
    echo "${manuscript_num},${pipeline_dirname},,no,no,0,0,,,,,,,,," >> "${OUT_CSV}"
    continue
  fi
  for s in "${siblings[@]}"; do
    sibling_basename=$(basename "$s")
    if ! probe_dir "$manuscript_num" "$pipeline_dirname" "$sibling_basename"; then
      log "  WARNING: probe_dir failed for ${sibling_basename} (exit $?), writing partial row"
      echo "${manuscript_num},${pipeline_dirname},${sibling_basename},,,0,0,,,,,,,,," >> "${OUT_CSV}"
    fi
  done
done

log ""
log "Inventory written to: ${OUT_CSV}"
log "Log written to:       ${OUT_LOG}"
log ""
log "=== Summary ==="
sites_with_chm=$(awk -F, 'NR>1 && $5=="yes" {print $1}' "${OUT_CSV}" | sort -u | wc -l | tr -d ' ')
total_chm_dirs=$(awk -F, 'NR>1 && $5=="yes"' "${OUT_CSV}" | wc -l | tr -d ' ')
sites_multi_year=$(awk -F, 'NR>1 && $5=="yes" {c[$1]++} END{n=0; for(k in c) if(c[k]>=2) n++; print n}' "${OUT_CSV}")
log "  Sites with at least one CHM directory found: ${sites_with_chm} / 19"
log "  Total CHM directories found (across all sites): ${total_chm_dirs}"
log "  Sites with >=2 CHM directories (multi-year candidates): ${sites_multi_year}"
log ""
log "Outputs:"
log "  ${OUT_CSV}"
log "  ${OUT_LOG}"
