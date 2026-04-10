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

## All-Source Migration (2026-04-10)

- [x] All 11 source deployments migrated to current main (d74cf93) with zarr v3 fixes, pixi.toml, updated scripts
- [x] All `aggregate_broad.yml` updated: `ome-zarr` (zarr v2) → `zarr>=3`
- [x] All profiles standardized: 15 cores, vram_mb=80000, mem_mb=300000
- [x] All diameters standardized: nucleus=23, cytosol=57 (source04 was 25/69, changed for consistency)
- [x] Source02 samples.json: removed subset filters (Metadata_PlateType, Metadata_InChIKey) — now processes full source_2
- [x] New zarr v3 conda env created (hash d07e64204b55eb5dcc9dfdcb6fbe5b55)

## Per-Source Pipeline Status

- [ ] **source02**: Migrated. Previously complete but needs aggregation rerun with zarr v3 for consistency.
- [ ] **source03**: 98.6% extracted (4128/4186). Migrated, not running — needs restart for last 58 batches, then aggregation.
- [ ] **source04**: Migrated. **Needs full rerun** — diameters changed from 25/69 to 23/57, all prior extraction/segmentation is invalid.
- [ ] **source05**: 80.5% extracted (2896/3595). Migrated, not running.
- [ ] **source06**: 100% extracted (3207). Aggregation pipeline running as of 2026-04-10 ~22:00. Currently in create_broad_structure phase.
- [ ] **source07**: 53.9% extracted (953/1769). Migrated, not running.
- [ ] **source08**: Migrated. Previously complete but needs aggregation rerun with zarr v3 for consistency.
- [ ] **source09**: 64.7% extracted (1717/2652). Migrated, not running.
- [ ] **source10**: Migrated. Previously complete but needs aggregation rerun with zarr v3 for consistency.
- [ ] **source11**: 73.7% extracted (1779/2415). Migrated, not running.
- [ ] **source13**: 100% extracted (286), 93 segmentations. Migrated, not running — needs pipeline run.
- [ ] Kick off remaining sources on separate nodes

## Source 04 Completed Items (prior run — now invalidated by diameter change)

- [x] Diagnosed why source_4 only had 244/708 batches: May 2025 re-run was killed mid-execution
- [x] Installed snakemake 8 via pixi, created pixi.toml, updated profile (15 cores, 80GB VRAM, 300GB RAM)
- [x] Pipeline completed successfully (2026-04-08 19:22) — all 277 plates rechunked
- [x] ~~Verify rechunked compressed output integrity~~ — moot, full rerun needed (diameters changed 25/69 → 23/57)
- [ ] Clean up old source04 results before rerun (extraction, segmentation, aggregated dirs contain data from old diameters)
- [ ] Clean up leftover temp mmap dirs in source04 `results/tmp/` (40 dirs from May 2025)

## Performance

- [x] Profiled bottleneck: HDF5 extraction over NFS (~16 min/batch) dominates; GPU cellpose takes only ~2 min
- [x] Raised profile limits: vram_mb 15k→70k, mem_mb 100k→240k
- [x] Discovered `mem_mb` resource limit doesn't cap snakemake concurrency — `cores` is the reliable lever
- [x] Set `cores: 12` as hard cap (12 parallel extract jobs confirmed working)

## Downstream

- [x] ~~Replicate zarr v3 fixes to other source pipelines~~ — all 11 sources migrated to current main (2026-04-10)
- [ ] Verify rechunk_broad Zarr v3 output is correct (check source06 output once aggregation finishes)
- [ ] Final aggregate step for all sources still pending
- [ ] Rerun aggregation for sources 02, 08, 10 with zarr v3 for output consistency

## Config Notes (all sources, NOT committed — local overrides in each deployment)
- `samples.json`: `{"samples": [{"Metadata_Source": "source_N"}]}` per source
- `sparcspy.yml`: nucleus diameter=23, cytosol diameter=57 (all sources standardized)
- `profile/config.yaml`: 15 cores, vram_mb=80000, mem_mb=300000 (all sources standardized)
- `aggregate_broad.yml`: zarr>=3 (all sources, replaces ome-zarr/zarr v2)

## Known Issues
- Extract rule outputs are marked `temp()` — extraction/segmentation H5 files get deleted after downstream rules consume them
- Snakemake warns about `directory` flag used in inputs for `aggregate_broad_batch` and `rechunk_broad` rules (non-fatal)
