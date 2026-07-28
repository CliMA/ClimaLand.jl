#!/bin/bash
# Summarize model-vs-MODIS LAI skill across one or more battery runs into a
# markdown table, for pasting into PR comments. Scans each given run directory
# for */out/lai_metrics.txt (written by single_site.jl) and prints one row per
# site, plus a mean-RMSE footer per run directory.
#
# Usage:
#   summarize_lai.sh <battery_rundir> [<battery_rundir2> ...]
# Example (compare with-beta vs beta-free):
#   summarize_lai.sh .../battery_6834627.desched1 .../battery_6834950.desched1
set -u

if [[ $# -lt 1 ]]; then
    echo "usage: $0 <battery_rundir> [<battery_rundir2> ...]" >&2
    exit 1
fi

echo "| run | site | beta | online_f0 | f0 | LAI_RMSE | LAI_BIAS | n |"
echo "|-----|------|------|-----------|----|----------|----------|---|"
for rundir in "$@"; do
    tag=$(basename "$rundir")
    sum=0
    cnt=0
    while IFS= read -r f; do
        site=$(awk '/^site /{print $2}' "$f")
        beta=$(awk '/^beta_in_A0 /{print $2}' "$f")
        onl=$(awk '/^online_f0 /{print $2}' "$f")
        f0=$(awk '/^f0 /{print $2}' "$f")
        rmse=$(awk '/^LAI_RMSE /{print $2}' "$f")
        bias=$(awk '/^LAI_BIAS /{print $2}' "$f")
        n=$(awk '/^n /{print $2}' "$f")
        printf "| %s | %s | %s | %s | %s | %.3f | %+.3f | %s |\n" \
            "$tag" "$site" "${beta:-?}" "${onl:-?}" "${f0:-?}" "$rmse" "$bias" "${n:-?}"
        sum=$(awk -v s="$sum" -v r="$rmse" 'BEGIN{print s + r}')
        cnt=$((cnt + 1))
    done < <(find "$rundir" -name lai_metrics.txt 2>/dev/null | sort)
    if [[ $cnt -gt 0 ]]; then
        mean=$(awk -v s="$sum" -v c="$cnt" 'BEGIN{printf "%.3f", s / c}')
        printf "| **%s** | **mean (%d sites)** | | | | **%s** | | |\n" \
            "$tag" "$cnt" "$mean"
    fi
done
