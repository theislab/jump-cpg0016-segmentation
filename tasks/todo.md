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
- [x] All diameters standardized: nucleus=23, cytosol=57 (all sources, including source04)
- [x] All sources pulled to main (2edee45) and pixi installed (2026-04-11)
- [x] Source02 samples.json: removed subset filters (Metadata_PlateType, Metadata_InChIKey) — now processes full source_2
- [x] New zarr v3 conda env created (hash d07e64204b55eb5dcc9dfdcb6fbe5b55)

## Per-Source Pipeline Status

**CRITICAL (discovered 2026-04-11): ALL 2024 extractions used wrong diameters (21/53 instead of 23/57). Every source needs full rerun except source04 and source13.**

- [ ] **source02**: Needs full rerun. 1960 batches. Old 2024 extraction invalid (21/53). No extraction dir remains.
- [ ] **source03**: Full rerun started on separate node (2026-04-11). 4186 batches. Config fixed to 23/57.
- [x] **source04**: **COMPLETE and ready for Broad transfer.** Extracted Apr 8 2026 with correct 23/57. 708 batches, 277 plates rechunked. aggregate.txt present.
- [ ] **source05**: Needs full rerun. 3595 batches. Old 2024 extraction invalid (21/53). Old results need cleaning.
- [ ] **source06**: Needs full rerun. 3206 batches. Old 2024 extraction invalid (21/53). Old aggregation was on invalid data.
- [ ] **source07**: Full rerun started on hpc-build01 (PID 4049212, 2026-04-11). Old results cleaned. 1769 batches.
- [ ] **source08**: Needs full rerun. 2986 batches. Old 2024 extraction invalid (21/53). No extraction dir remains.
- [ ] **source09**: Needs full rerun. 2652 batches. Old 2024 extraction invalid (21/53). Old results need cleaning.
- [ ] **source10**: Needs full rerun. 2038 batches. Old 2024 extraction invalid (21/53). No extraction dir remains.
- [ ] **source11**: Needs full rerun. 2415 batches. Old 2024 extraction invalid (21/53). Old results need cleaning.
- [ ] **source13**: Running on gpusrv65 with correct 23/57. ~286 batches extracted as of 2026-04-11.
- [ ] Kick off remaining sources (02, 05, 06, 08, 09, 10, 11) on separate nodes — clean old results first
- [ ] Transfer source04 to Broad (first source ready)

## Source 04 — COMPLETE (Apr 8 2026, correct 23/57 diameters)

- [x] Diagnosed why source_4 only had 244/708 batches: May 2025 re-run was killed mid-execution
- [x] Installed snakemake 8 via pixi, created pixi.toml, updated profile (15 cores, 80GB VRAM, 300GB RAM)
- [x] Pipeline completed successfully (2026-04-08 19:22) — all 277 plates rechunked with correct 23/57 diameters
- [x] Verified: extraction ran Apr 8 2026 (after diameter fix), aggregate.txt present, rechunk checkpoint complete
- [ ] Clean up leftover temp mmap dirs in source04 `results/tmp/` (40 dirs from May 2025)

## Performance

- [x] Profiled bottleneck: HDF5 extraction over NFS (~16 min/batch) dominates; GPU cellpose takes only ~2 min
- [x] Raised profile limits: vram_mb 15k→70k, mem_mb 100k→240k
- [x] Discovered `mem_mb` resource limit doesn't cap snakemake concurrency — `cores` is the reliable lever
- [x] Set `cores: 12` as hard cap (12 parallel extract jobs confirmed working)

## Downstream

- [x] ~~Replicate zarr v3 fixes to other source pipelines~~ — all 11 sources migrated to current main (2026-04-10)
- [ ] Verify rechunk_broad Zarr v3 output is correct (check source04 output — first complete source)
- [ ] Transfer source04 compressed zarr output to Broad
- [ ] All sources need full pipeline completion (extract → aggregate → rechunk) with correct 23/57 diameters

## Config Notes (all sources, NOT committed — local overrides in each deployment)
- `samples.json`: `{"samples": [{"Metadata_Source": "source_N"}]}` per source
- `sparcspy.yml`: nucleus diameter=23, cytosol diameter=57 (all sources standardized)
- `profile/config.yaml`: 15 cores, vram_mb=80000, mem_mb=300000 (all sources standardized)
- `aggregate_broad.yml`: zarr>=3 (all sources, replaces ome-zarr/zarr v2)

## Known Issues
- Extract rule outputs are marked `temp()` — extraction/segmentation H5 files get deleted after downstream rules consume them
- Snakemake warns about `directory` flag used in inputs for `aggregate_broad_batch` and `rechunk_broad` rules (non-fatal)
