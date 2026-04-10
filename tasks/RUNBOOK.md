# JUMP CPG0016 Segmentation — Pipeline Runbook

Self-contained instructions for running the segmentation pipeline on GPU cluster nodes.
Last updated: 2026-04-10.

---

## Quick Reference

| Source | Batches | Extracted | Status | What To Do |
|--------|---------|-----------|--------|------------|
| 02 | TBD (full source, no filters) | previously complete | needs full rerun | old run used subset filters; rerun from scratch for full source_2 |
| 03 | 4186 | 4128 (98.6%) | not running | finish last 58 batches + aggregation |
| 04 | 708 | previously complete | **needs full rerun** | diameters changed (was 25/69, now 23/57); old results invalid |
| 05 | 3595 | 2896 (80.5%) | not running | continue extraction + aggregation |
| 06 | 3207 | 3207 (100%) | **aggregation running** | started 2026-04-10 ~22:00; monitor |
| 07 | 1769 | 953 (53.9%) | not running | continue extraction + aggregation |
| 08 | TBD | previously complete | needs aggregation rerun | rerun aggregation with zarr v3 for consistency |
| 09 | 2652 | 1717 (64.7%) | not running | continue extraction + aggregation |
| 10 | TBD | previously complete | needs aggregation rerun | rerun aggregation with zarr v3 for consistency |
| 11 | 2415 | 1779 (73.7%) | not running | continue extraction + aggregation |
| 13 | 286 | 286 (100%) | not running | run pipeline (93 segmentations exist, rest will be filled) |

### Priority Order

1. **source06** — already running, just monitor
2. **source13** — 286 batches, smallest source, quick win
3. **source03** — only 58 batches left, then aggregation
4. **source08/10** — aggregation-only reruns (extraction already done)
5. **source02** — full rerun (old run used subset filters)
6. **source05/09/11** — partially extracted, continue
7. **source04/07** — largest remaining work

---

## How To Run a Source

### 1. Get on a GPU node

You need a node with:
- 1x GPU (any NVIDIA, A100 preferred but not required — GPU is NOT the bottleneck)
- 15+ CPU cores
- 300+ GB RAM
- NFS access to `/ictstr01/groups/ml01/`

### 2. Navigate to the source directory

```bash
cd /ictstr01/groups/ml01/projects/2023_ttreis_segment_JUMP/snakemake/final_source{NN}/snakemake
```

Replace `{NN}` with the source number (02, 03, 04, ..., 13).

### 3. Verify the deployment is on current main

```bash
git log --oneline -1
# Should show: d74cf93 Merge pull request #10 from theislab/add-tasks-and-lessons
```

If it doesn't, the source hasn't been migrated. See "Migration" section below.

### 4. Check for stale locks

```bash
ls .snakemake/locks/
```

If there are `.lock` files:
```bash
# First check no snakemake is actually running
pgrep -af "snakemake.*final_source{NN}"

# If nothing running, clear locks
find .snakemake/locks/ -name "*.lock" -delete
```

### 5. Dry run (optional but recommended for first time)

```bash
pixi run dryrun
```

This builds the DAG (~30 min for large sources on NFS) and shows what jobs will run. No actual work is done.

### 6. Start the pipeline

```bash
pixi run run
```

This runs in the foreground. For background execution:
```bash
nohup pixi run run > ../log/run_$(date +%Y%m%d_%H%M%S).log 2>&1 &
echo $! > ../log/snakemake.pid
```

### 7. Monitor progress

**Do NOT trust log files** — NFS write cache causes 30-60+ minute lag.

Count output files on disk instead:
```bash
# Extractions done
ls results/extraction/source_{N}/ | wc -l

# Segmentations done
ls results/segmentation/source_{N}/ | wc -l

# Aggregate batch checkpoints done
ls results/checkpoints/aggregate_broad_batch/source_{N}/ 2>/dev/null | wc -l

# Check if pipeline process is alive
pgrep -af "snakemake.*final_source{NN}"
```

Or use the progress script from the repo root:
```bash
bash /ictstr01/groups/ml01/workspace/ttreis/projects/jump-cpg0016-segmentation/scripts/check_progress.sh
```

---

## Special Cases

### Sources needing FULL rerun (02, 04)

For source04 (diameter change) and source02 (old subset filters), old results must be cleaned first:

```bash
cd /ictstr01/groups/ml01/projects/2023_ttreis_segment_JUMP/snakemake/final_source{NN}/snakemake

# Check what exists
ls results/

# Remove old extraction/segmentation/aggregation results
# WARNING: This deletes data. Only do this for sources that need full rerun.
rm -rf results/extraction/ results/segmentation/ results/aggregated/ results/checkpoints/
rm -rf results/tmp/  # old temp mmap dirs

# Then run normally
pixi run run
```

### Sources needing aggregation-only rerun (08, 10)

These have valid extraction+segmentation but need aggregation redone with zarr v3:

```bash
cd /ictstr01/groups/ml01/projects/2023_ttreis_segment_JUMP/snakemake/final_source{NN}/snakemake

# Remove old aggregation outputs only
rm -rf results/aggregated/ results/checkpoints/

# Snakemake will see extraction+segmentation exist and only run aggregation rules
pixi run run
```

### Source06 (currently running)

Already running as of 2026-04-10 ~22:00. Pipeline steps:
1. `create_broad_structure` — single job, iterates all 3207 batches to create dir layout
2. `aggregate_broad_batch` × 3207 — per-batch zarr aggregation (15 parallel)
3. `rechunk_broad` — compress all plates into sharded zarr v3
4. `aggregate` — write final checkpoint to `results/aggregated/aggregate.txt`

Check if still running:
```bash
pgrep -af "snakemake.*final_source06"
```

---

## Standardized Config (all sources)

Every source deployment has these local overrides (not committed to git):

**`config/samples.json`**:
```json
{"samples": [{"Metadata_Source": "source_N"}]}
```

**`config/sparcspy.yml`** (diameters):
- nucleus: 23
- cytosol: 57

**`profile/config.yaml`**:
```yaml
cores: 15
use-conda: True
keep-incomplete: True
rerun-incomplete: True
resources:
  vram_mb: 80000
  mem_mb: 300000
```

**`envs/aggregate_broad.yml`**: uses `zarr>=3` (NOT `ome-zarr`)

---

## Troubleshooting

### Pipeline hangs after "Shutting down"
Snakemake gets stuck in multiprocessing cleanup on NFS. Kill it:
```bash
kill -9 $(pgrep -f "snakemake.*final_source{NN}")
# Also kill orphaned child processes
ps aux | grep "scripts/extract.py" | grep -v grep | awk '{print $2}' | xargs kill -9 2>/dev/null
# Clear locks before restarting
find .snakemake/locks/ -name "*.lock" -delete
```

### Lock files blocking startup
```bash
find .snakemake/locks/ -name "*.lock" -delete
```
Don't use `pixi run unlock` — it rebuilds the full DAG (~30 min) just to unlock.

### DAG build takes forever
Normal. 30-35 min for ~6000+ jobs on NFS. The pipeline isn't hung, it's just building.

### "TiffFileError: not a TIFF file: header=b''"
Transient S3 empty response. The rule has `retries: 5`. If all fail, they get fresh retries on restart.

### conda env creation fails
The shared conda prefix is at:
```
/ictstr01/groups/ml01/projects/2023_hackathon23_subcellular_spatial_niklas.schmacke/jump/envs/
```
Key env hashes:
- Download (boto3): `507f1bacc8e964624ea0c0fc2f042010_`
- Extract (SPARCSpy): `79ff5f567d72279473b1703565acfd26_`
- Aggregate broad (zarr v3): `d07e64204b55eb5dcc9dfdcb6fbe5b55_`
- Rechunk (zarr v3 + dask): `67cb26037c82092ba47126f5b5e03793_`

### GPU is idle most of the time
Normal. The bottleneck is HDF5 extraction over NFS (~16 min/batch), not GPU cellpose (~2 min/batch). GPU utilization ~15%.

---

## Pipeline Architecture

```
download_stacks → extract → create_broad_structure → aggregate_broad_batch (×N) → rechunk_broad → aggregate
```

- **download_stacks**: Downloads .npy stacks from S3 (cellpainting-gallery)
- **extract**: Runs SPARCSpy + Cellpose segmentation, extracts single cells to H5
- **create_broad_structure**: Creates directory layout for Broad-format output
- **aggregate_broad_batch**: Per-batch: reads H5, writes zarr groups (label_image, single_cell_index, single_cell_data)
- **rechunk_broad**: Recompresses all zarr plates with sharded zarr v3 (zstd, per-cell chunks)
- **aggregate**: Writes final checkpoint file

Extraction outputs are marked `temp()` — H5 files are deleted after aggregation consumes them.

---

## Running Multiple Sources in Parallel

Each source needs its own GPU node (or at least its own GPU + 15 cores). To run N sources simultaneously:

1. Get N GPU nodes
2. On each node, `cd` to a different `final_source{NN}/snakemake/`
3. Run `pixi run run` on each

Sources are completely independent — they don't share state. The only shared resource is the conda env prefix directory, which is read-only after env creation.
