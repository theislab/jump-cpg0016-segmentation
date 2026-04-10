#!/usr/bin/env bash
# Check segmentation pipeline progress for all sources.
# Counts output files on disk (not log messages) to avoid NFS cache lag.

SNAKEMAKE_BASE=/ictstr01/groups/ml01/projects/2023_ttreis_segment_JUMP/snakemake
PYTHON=/home/icb/tim.treis/miniforge3/bin/python

printf "%-10s %-8s %-8s %-8s %-8s %-12s\n" "Source" "Extracted" "Total" "Pct%" "Images" "Broad_agg"
printf "%-10s %-8s %-8s %-8s %-8s %-12s\n" "------" "---------" "-----" "----" "------" "---------"

for src in 02 03 04 05 06 07 08 09 10 11 13; do
  src_num=${src#0}  # strip leading zero for directory names (source_3 not source_03)
  BASE=$SNAKEMAKE_BASE/final_source${src}/snakemake

  [ -d "$BASE" ] || continue

  extracted=$(ls "$BASE/results/extraction/source_${src_num}/" 2>/dev/null | wc -l)
  images=$(ls "$BASE/results/images/" 2>/dev/null | wc -l)
  broad=$(ls "$BASE/results/aggregated/broad/" 2>/dev/null | wc -l)

  total=$($PYTHON -c \
    "import pandas as pd; df=pd.read_parquet('$BASE/config/selected_metadata.parquet'); print(df['snakemake_batch'].nunique())" \
    2>/dev/null || echo "?")

  if [ "$total" = "?" ] || [ "$total" = "0" ]; then
    pct="?"
  else
    pct=$(python3 -c "print(f'{$extracted/$total*100:.0f}%')" 2>/dev/null || echo "?")
  fi

  final=""
  [ -f "$BASE/results/aggregated/aggregate.txt" ] && final=" [DONE]"

  printf "%-10s %-8s %-8s %-8s %-8s %-12s%s\n" \
    "source${src}" "$extracted" "$total" "$pct" "$images" "$broad" "$final"
done

echo ""
echo "Running snakemake processes:"
ps aux | grep "snakemake.*profile" | grep -v grep | awk '{print "  PID "$2" ("$11")"}' || echo "  none"
