# JUMP CPG0016 Segmentation — Task Tracker

## Pipeline Fixes (snakemake repo: theislab/jump-cpg0016-segmentation)

- [x] Fix `download_stacks` rule: replaced `conda:` + `script:` with `shell:` using direct conda env Python path (`507f1bacc8e964624ea0c0fc2f042010_`)
  - Files: `snakemake/rules/download.smk`, `snakemake/scripts/download_stack.py`
- [x] Fix `extract` rule: same pattern — replaced `conda:` + `script:` with `shell:` using direct conda env Python path (`79ff5f567d72279473b1703565acfd26_`)
  - Files: `snakemake/rules/extract.smk`, `snakemake/scripts/extract.py`
- [x] Add argparse CLI to `download_stack.py` so it can be called from `shell:` (previously only had magic `snakemake` object support)
- [x] Add argparse CLI to `extract.py` (same reason)
- [x] Add rechunk_broad rule for sharded Zarr v3 compression (separate prior commit)
- [x] All above pushed to PR #9 (`feat/rechunk-snakemake-rule`): https://github.com/theislab/jump-cpg0016-segmentation/pull/9
- [x] PR #9 merged to main
- [x] Fixed zarr v3 incompatibility in `aggregate_broad_batch.py`: removed context manager on `zarr.open()`, changed `.array()` → `.create_array()` (PR #11)
- [x] Fixed zarr v3 incompatibility in `rechunk_broad.py`: replaced dask `to_zarr()` with native zarr v3 `create_array()` using `shards`/`chunks`/`compressors`; initialized store before workers (PR #11)
- [x] Added pixi.toml for snakemake 8 (PR #11)

## Per-Source Pipeline Status

- [x] sources 02, 08, 10: Fully complete (final aggregate exists).
- [x] **source04**: Completed 2026-04-08. All 708 batches extracted, 277 plates rechunked, aggregate checkpoint written.
- [ ] **source03**: 98.6% done (4128/4186 batches extracted). Not currently running — needs restart for last 58 batches, then aggregation.
- [ ] **source06**: 100% extraction done (3207/3206), 2 broad batches aggregated — needs aggregation pipeline run.
- [ ] **source13**: 100% extraction done (286/286 batches) — needs aggregation pipeline run.
- [ ] **source05**: 80.5% done (2896/3595 batches). Not running.
- [ ] **source09**: 64.7% done (1717/2652 batches). Not running.
- [ ] **source11**: 73.7% done (1779/2415 batches). Not running.
- [ ] **source07**: 53.9% done (953/1769 batches). Not running.
- [ ] Kick off remaining sources on separate nodes

## Source 04 Completed Items

- [x] Diagnosed why source_4 only had 244/708 batches: May 2025 re-run was killed mid-execution
- [x] Installed snakemake 8 via pixi, created pixi.toml, updated profile (15 cores, 80GB VRAM, 300GB RAM)
- [x] Pipeline completed successfully (2026-04-08 19:22) — all 277 plates rechunked
- [ ] Verify rechunked compressed output integrity (spot-check a few plates)
- [ ] Clean up leftover temp mmap dirs in source04 `results/tmp/` (40 dirs from May 2025)
- [ ] Clean up duplicate log files from failed launch attempts

## Performance

- [x] Profiled bottleneck: HDF5 extraction over NFS (~16 min/batch) dominates; GPU cellpose takes only ~2 min
- [x] Raised profile limits: vram_mb 15k→70k, mem_mb 100k→240k
- [x] Discovered `mem_mb` resource limit doesn't cap snakemake concurrency — `cores` is the reliable lever
- [x] Set `cores: 12` as hard cap (12 parallel extract jobs confirmed working)

## Downstream

- [ ] Replicate zarr v3 fixes to other source pipelines if they use older scripts
- [ ] Verify rechunk_broad Zarr v3 output is correct before running on all sources
- [ ] Final aggregate step for all sources still pending

## Config Notes (source_4-specific, NOT committed to main)
- `samples.json`: `{"samples": [{"Metadata_Source": "source_4"}]}` (708 batches)
- `sparcspy.yml`: nucleus diameter=25, cytosol diameter=69 (differs from repo default of 21/53)
- `profile/config.yaml`: 15 cores, vram_mb=80000, mem_mb=300000

## Known Issues
- Extract rule outputs are marked `temp()` — extraction/segmentation H5 files get deleted after downstream rules consume them
- Snakemake warns about `directory` flag used in inputs for `aggregate_broad_batch` and `rechunk_broad` rules (non-fatal)
