# Refactor: Aggregation to per-Well Grouping (Option B)

Status: rule + script implemented and CLI-validated on one plate, 2026-05-20.

## Motivation

Switch zarr output from per-FOV groups (Option A: ~24k groups/plate → ~60M total objects across cpg0016) to per-well groups with FOV concatenation (Option B: ~3.5k groups/plate → ~8.6M total). Broad team (Erin, shntnu) prefers B because well is the unit of perturbation and the S3 LIST/sync cost on tens of millions of objects is prohibitive. While we rewrite, also fix the provenance gap — current `single_cell_index` ([index, cellid]) is meaningless once flattened from group hierarchy.

## Strategy

Do **not** rewrite the existing snakemake aggregation rules. The pipeline keeps producing Option A. Add a **new** rule `reaggregate_to_well` that reads one Option A `{plate}.zarr` and writes one Option B `{plate}.zarr` in a sibling tree. The rule is the unit of fan-out: one job per plate, snakemake handles parallelism, retries, restart. This:

- Decouples the schema migration from the running pipelines (6 sources currently in flight).
- Treats the 5 completed sources (02, 04, 10, 11, 13) and the 6 in-flight ones uniformly — same rule applies after each plate's Option A is committed.
- Loses zero provenance. `image_id` and `label_id` are recoverable from the source zarr's group names and `single_cell_index[:, 1]` respectively; `selected_metadata.parquet` is not needed at conversion time.

## Final schema

```
{plate}.zarr/                                [group]
  attrs:
    source        e.g. "source_4"            (parsed from any group name)
    batch         e.g. "Batch1"              (parsed from any group name)
    plate         e.g. "BR00121439"          (parsed from any group name)
    fov_height    e.g. 1080                  (from src label_image.shape[-2])
    fov_width     e.g. 1080                  (from src label_image.shape[-1])
    n_wells       e.g. 384                   (from grouping)
    created_at    ISO-8601 UTC

  channel_mapping.json                       [sibling file, copied verbatim from Option A]

  {WELL}/                                    [group, e.g. "A01"]
    attrs:
      well     "A01"
      n_fovs   e.g. 9
      n_cells  e.g. 3741
      fov_lut  list[str]                    # one full image_id per FOV slot:
                                            # fov_lut[fov_idx] -> original
                                            # "source__batch__plate__well__fov"
                                            # Stored as a JSON attr (not a zarr
                                            # array) because zarr v3 has no
                                            # stable spec for fixed-length
                                            # unicode; attrs are always portable.

    label_image        uint32   (n_fovs, 2, H, W)
                       chunks=(1, 1, H, W)
                       shards=(n_fovs, 2, H, W)
                       compressor=BloscCodec(zstd, clevel=3, shuffle)
                       # (fov_idx, label_id) is the unique key for any cell.

    single_cell_data   float16  (N_cells, 7, 150, 150)
                       chunks=(1, 7, 150, 150)
                       shards=(min(4096, N_cells), 7, 150, 150)
                       compressor=BloscCodec(zstd, clevel=3, shuffle)

    single_cell_index  uint32   (N_cells, 2)
                       cols = [fov_idx, label_id]
                       # fov_idx  → axis-0 of label_image / fov_lut
                       # label_id → label present in label_image[fov_idx]
```

**Round-trip guarantee:** for any cell row `k`, `mask = label_image[fov_idx, channel] == label_id` is non-empty, where `(fov_idx, label_id) = single_cell_index[k]`.

### Locked schema decisions

| Decision | Choice |
|---|---|
| Output location | sibling tree `aggregated_v2/` — Option A stays untouched |
| `fov_lut` storage | well-group attr (JSON list of full `image_id` strings) — avoids zarr v3 unicode-spec instability |
| `single_cell_index` columns | 2 cols `[fov_idx, label_id]` — `snakemake_batch_idx` and `src_row` dropped for cleaner Broad delivery |
| `channel_mapping.json` | sibling file only (not duplicated into attrs) — single source of truth |
| `pipeline_version` attr | dropped — `created_at` is sufficient provenance |
| `selected_metadata.parquet` at conversion time | not used — everything derivable from source zarr |
| Label-id relabeling on concat | none — labels stay FOV-local |
| Variable `n_fovs` per well | yes, sized to actual count per well (no padding) |

## Path mapping

```
src: results/aggregated/broad/cellpainting-gallery/cpg0016-jump/{source}/workspace/segmentation/{model_dir}/objects/{batch}/{plate}/{plate}.zarr
dst: results/aggregated_v2/cellpainting-gallery/cpg0016-jump/{source}/workspace/segmentation/{model_dir}/objects/{batch}/{plate}/{plate}.zarr
```

The pipeline pivoted away from `broad_compressed/` entirely. `rule rechunk_broad` is disabled (code preserved). The reaggregator reads from `broad/` directly, which is verified byte-equal to `broad_compressed/` for all done sources — including source04, where the `cellpose/` vs `cellpose_202404/` directory rename was cosmetic-only (same underlying segmentation data).

Wall-clock impact of reading the unsharded `broad/` instead of sharded `broad_compressed/`: none observed on the one-plate test (171 s either way). Reaggregation is sequential per FOV, so sharding gave no measurable benefit.

`{model_dir}` is a wildcard for forward-compatibility (any future model-version subdirectory just flows through verbatim). For all currently-existing sources, `{model_dir}` resolves to `cellpose`.

## Snakemake rule

Added to `rules/aggregate.smk`:

```python
rule reaggregate_to_well:
    input:
        plate_zarr=directory(".../broad/.../{source}/.../{model_dir}/objects/{batch}/{plate}/{plate}.zarr"),
        # All per-batch aggregations for the source must be committed before any
        # plate is reaggregated (a plate's FOVs may have been written by multiple
        # snakemake batches into the same broad/{plate}.zarr).
        aggregate_batch_checkpoints=get_aggregate_broad_batch_checkpoints,
    output:
        checkpoint="results/checkpoints/reaggregate_to_well/{source}/{model_dir}/{batch}/{plate}.ckpt",
        # Tracks the plate-level *parent* dir, not the .zarr itself, so snakemake's
        # .snakemake_timestamp marker file lives next to the zarr rather than inside it
        # (zarr readers warn on unrecognized entries at the zarr root).
        plate_v2_dir=directory(".../aggregated_v2/.../{source}/.../{model_dir}/objects/{batch}/{plate}"),
    log: "log/reaggregate_to_well/{source}/{model_dir}/{batch}/{plate}.log"
    conda: "../envs/aggregate_broad.yml"
    resources:
        mem_mb=4000,
        runtime=120,
    retries: 2
    script: "../scripts/reaggregate_to_well.py"
```

`rule rechunk_broad` is commented out (code preserved in `rules/aggregate.smk`). The `broad_compressed/` tree it produced is no longer a pipeline output. Existing `broad_compressed/` data on disk for done sources can be deleted whenever convenient — see "Cleanup" below.

### Barrier rule: `aggregate_broad_done`

Per-plate `reaggregate_to_well` jobs need every per-batch `aggregate_broad_batch` to be committed first (a plate's FOVs can span multiple snakemake batches). Wiring each plate to all N batch checkpoints would mean ~N × (plates/source) input edges (~44k stats for source04) — slow to plan and fragile. Added a one-line barrier rule:

```python
rule aggregate_broad_done:
    input:
        batch_checkpoints=get_aggregate_broad_batch_checkpoints,
    output:
        checkpoint="results/checkpoints/aggregate_broad_done/{source}.ckpt",
    run:
        ...write a marker file...
```

`reaggregate_to_well` then carries a single source-level dependency (`aggregate_broad_done/{source}.ckpt`) and one per-plate directory check. DAG construction collapses from ~44k stats to a few hundred.

### Final target wiring

`rule aggregate` in the Snakefile now expands `ALL_REAGGREGATE_CHECKPOINTS` (every `(source, batch, plate)` tuple from `selected_metadata.parquet`, with `model_dir` hardcoded to `cellpose`). Output `aggregated/aggregate.txt` summarises plate counts per source rather than dumping every path (would hit `ARG_MAX`).

The rule is not wired into `rule all` yet — invoke it by requesting a specific checkpoint, or add a new top-level target after the one-plate test is approved.

The script (`scripts/reaggregate_to_well.py`) also exposes an argparse CLI for direct standalone use (matches the dual-entry pattern of `extract.py`):

```bash
reaggregate_to_well.py \
  --src-plate-zarr  <path to Option A {plate}.zarr> \
  --dst-plate-zarr  <path to Option B {plate}.zarr> \
  [--validate-sample 5] \
  [--wells A01,A02]                       # optional, for testing only
```

## Converter algorithm

Per plate, in order:

1. Open source zarr read-only. Enumerate FOV subgroups. Parse each name as `source__batch__plate__well__fov`. Group by well; numerically sort each well's FOVs by integer FOV (NOT lex — `"10"` must sort after `"2"`).
2. Reject if `{plate}.zarr` already exists at destination. Wipe any leftover `{plate}.zarr.tmp` from a prior crash.
3. Create `{plate}.zarr.tmp` as a zarr v3 group. Set plate-level attrs by deriving from the first group's image_id + first label_image shape.
4. For each well (serial):
   - Compute `n_fovs`, `n_cells = sum(per-FOV cell counts)`.
   - Create `label_image`, `fov_lut`, `single_cell_data`, `single_cell_index` with the shapes/chunks/shards/compressors above.
   - For `fov_idx, group_name in enumerate(sorted_fovs)`:
     - `label_image[fov_idx] := src[group_name].label_image`
     - Append `src[group_name].single_cell_data` at running cell offset.
     - Append `[fov_idx, src[group_name].single_cell_index[:, 1]]` at running cell offset. (Column 1 in Option A is `cellid` = the FOV-local `label_id`, verified by inspecting an extraction h5 in this session.)
     - `fov_lut[fov_idx] := group_name`
   - Set well group attrs.
5. Copy sibling `channel_mapping.json` verbatim into the v2 plate's parent dir.
6. **Validate before commit:** sample N random cells (default 5) across random wells; for each, read `(fov_idx, label_id)` from `single_cell_index`, fetch `mask = label_image[fov_idx, 1] == label_id`, assert `mask.any()`. Failure → exit non-zero, leave `.tmp/` for inspection.
7. Atomic `rename({plate}.zarr.tmp, {plate}.zarr)`. The plate is now committed.

## Atomicity & idempotency

- Write under `{plate}.zarr.tmp/`, single `rename` to `{plate}.zarr/` at the very end. The directory rename is atomic.
- Crash leaves `.tmp/` behind; next run wipes it before retrying.
- Refuse to start if `{plate}.zarr/` already exists — never overwrite a committed plate.
- No per-well checkpointing inside a plate. A plate is the unit of restart.

## Memory & parallelism

Per-well peak working set: ~9 FOVs × 2 × 1080² × 4 B (label_image) + ~9 × ~500 cells × 7 × 150² × 2 B (single_cell_data) ≈ **~220 MB**.

→ One well at a time in memory. Plate-level parallelism comes from running multiple plate jobs in parallel via SLURM array. No threading inside the script. If a benchmark shows a single plate needs inner parallelism, add a `--workers` flag for per-well `ProcessPoolExecutor` later.

## Implementation status

- [x] `scripts/reaggregate_to_well.py` (snakemake + CLI dual entry).
- [x] `rule reaggregate_to_well` added to `rules/aggregate.smk`.
- [x] CLI smoke test on 1 well of `source_4/2021_04_26_Batch1/BR00117035` — 15 s wall, round-trip OK.
- [x] CLI full-plate test on 18-well partial plate `BR00117035` — 159 s wall, 1.1 GB output, round-trip OK.
- [x] Snakemake invocation on the same plate — 171 s wall, all outputs landed, `.snakemake_timestamp` lives next to the zarr (parent-dir output), zarr opens warning-free.
- [ ] SLURM array on smallest done source after rule validation. Hand `sbatch` to user (per [[no-autonomous-slurm-submissions]]).
- [ ] Array-submit remaining done sources.

Scale extrapolation from the 18-well partial plate: ~9 s/well linearly → full 384-well plate ≈ 56 min, 23 GB. ~216 plates in source04 → ~3 h wall clock with a 50-job array.

## Cleanup (after full v2 reaggregation lands)

- `broad_compressed/` is now an orphan tree on disk for done sources (sources 02, 04, 10, 11, 13). Safe to delete — no rule reads from it. Batch as one SLURM cpu_p job per source per the 5.5h NFS rm note in `tasks/lessons.md`.
- `broad/` becomes deletable per source once `aggregated_v2/` has been produced and verified for every plate of that source. Until then, `broad/` is the only input to v2; do not delete prematurely.

## Out of scope (explicitly deferred)

- Wiring `reaggregate_to_well` into `rule all` / the default target. `rule aggregate` in the Snakefile still expects `rechunk_broad/{source}.ckpt`; for in-flight sources this would break a fresh `snakemake all`. Update to depend on per-plate `reaggregate_to_well` checkpoints when ready to scale.
- Wrapping `broad/` as `temp(directory(...))` in `rule create_broad_structure` so future runs auto-clean. One-line change, defer until after the rule has been validated across all done sources.
- Deleting old `aggregated/broad{,_compressed}/`. Manual + batched after each source's v2 lands.
- Updates to `tasks/todo.md` or pipeline-status memory.
