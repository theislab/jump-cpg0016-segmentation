# Source 04 Pipeline TODO

## Pipeline Run (2026-04-07)
- [x] Diagnosed why source_4 only had 244/708 batches: May 2025 re-run was killed mid-execution; original Feb 2024 run had completed all 708
- [x] Deleted `results_old/` (completed but outdated first run with wrong diameters)
- [x] Switched repo to `main` branch (reset to `origin/main`, PR #9 merged) while preserving source_4-specific config changes
- [x] Fixed snakemake — installed snakemake 8 via pixi (system snakemake 7 had broken pulp dependency)
- [x] Created `pixi.toml` for source_04 (mirroring source_03 setup)
- [x] Updated `profile/config.yaml` for current SLURM allocation: 15 cores, 80GB VRAM, 300GB RAM
- [x] Pipeline running: 1652 total steps (464 downloads + 476 extracts + 708 aggregate_broad_batch + create_broad_structure + rechunk_broad + aggregate)
- [x] Monitor pipeline to completion — all 1652 download/extract jobs done, 708 aggregate_broad_batch completed
- [x] Resume pipeline on new SLURM node
- [x] Fixed zarr v3 incompatibility in `aggregate_broad_batch.py`: removed context manager on `zarr.open()`, changed `.array()` → `.create_array()`
- [x] Fixed PATH ordering in `pixi.toml`: miniforge3/bin moved to end of PATH so conda env Python (3.14) takes precedence over base (3.12)
- [x] Fixed zarr v3 incompatibility in `rechunk_broad.py`: replaced dask `to_zarr()` with native zarr v3 `create_array()` using `shards`/`chunks`/`compressors`; removed `zarr_version=3` param; initialized store with `mode="w"` before spawning workers
- [x] Pipeline completed successfully (2026-04-08 19:22) — all 277 plates rechunked, aggregate checkpoint written
- [ ] Verify rechunked compressed output integrity (spot-check a few plates)
- [ ] Clean up leftover temp mmap dirs in `results/tmp/`
- [ ] Clean up duplicate log files in `../log/` from failed launch attempts (output_202604071139 through 202604071154)
- [ ] Replicate zarr v3 fixes to other source pipelines (source_02, source_03, source_08, source_13) if they use the same scripts

## Config Notes (source_4-specific, NOT committed to main)
- `samples.json`: `{"samples": [{"Metadata_Source": "source_4"}]}` (selects all 176,845 samples = 708 batches)
- `sparcspy.yml`: nucleus diameter=25, cytosol diameter=69 (differs from repo default of 21/53)
- `profile/config.yaml`: 15 cores, vram_mb=80000, mem_mb=300000

## Known Issues
- 12 batches were previously incomplete (64, 124, 125, 245, 364, 365, 454, 484, 485, 544, 574, 604) — being re-run with `--rerun-incomplete`
- Extract rule outputs are marked `temp()` — extraction/segmentation H5 files get deleted after downstream rules consume them
- Snakemake warns about `directory` flag used in inputs for `aggregate_broad_batch` and `rechunk_broad` rules (non-fatal)
