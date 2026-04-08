# Source 04 Pipeline Handover

## Status as of 2026-04-08 ~02:00
- **883/1652 jobs done (53%)**, 0 errors
- 689/708 images downloaded, 689 extraction dirs, 689 segmentation dirs
- ~19 remaining downloads, ~19 remaining extracts, 708 aggregate_broad_batch, 1 create_broad_structure, 1 rechunk_broad, 1 aggregate
- Pipeline will be killed when SLURM node dies (~2 min left)

## Folder Paths
- **Pipeline root:** `/ictstr01/groups/ml01/projects/2023_ttreis_segment_JUMP/snakemake/final_source04/`
- **Snakemake dir (run from here):** `/ictstr01/groups/ml01/projects/2023_ttreis_segment_JUMP/snakemake/final_source04/snakemake/`
- **Pixi env:** `/ictstr01/groups/ml01/projects/2023_ttreis_segment_JUMP/snakemake/final_source04/snakemake/.pixi/`
- **Pipeline logs:** `/ictstr01/groups/ml01/projects/2023_ttreis_segment_JUMP/snakemake/final_source04/log/`
- **Active log (this run):** `output_202604071156_stderr.log`
- **Extract logs:** `snakemake/log/extract/source_4/`
- **Config:** `snakemake/config/` (samples.json = source_4 only, sparcspy.yml = diam 25/69)
- **Profile:** `snakemake/profile/config.yaml` (15 cores, 80GB VRAM, 300GB RAM)
- **Results:** `snakemake/results/{images,extraction,segmentation,sparcspy,tmp}/`
- **Conda envs (shared, do not delete):** `/ictstr01/groups/ml01/projects/2023_hackathon23_subcellular_spatial_niklas.schmacke/jump/envs/`
- **Git repo:** `theislab/jump-cpg0016-segmentation`, branch `add-pixi-and-source04-config`, PR #11

## Resume Steps
1. Request new SLURM node (needs GPU, >=300GB RAM, >=15 cores):
   ```bash
   srun --partition=gpu_p --qos=gpu_normal --gres=gpu:1 --mem=300G --cpus-per-task=15 --time=1-00:00:00 --pty bash
   ```
2. cd to snakemake dir:
   ```bash
   cd /ictstr01/groups/ml01/projects/2023_ttreis_segment_JUMP/snakemake/final_source04/snakemake
   ```
3. Unlock (the kill will leave a stale lock):
   ```bash
   pixi run unlock
   ```
4. Adjust profile if new node has different resources:
   ```bash
   cat profile/config.yaml  # check cores/mem/vram values match your allocation
   ```
5. Run:
   ```bash
   timestamp=$(date +"%Y%m%d%H%M")
   pixi run run > "../log/output_${timestamp}_stdout.log" 2> "../log/output_${timestamp}_stderr.log" &
   ```
6. Verify it started:
   ```bash
   sleep 30 && tail -20 ../log/output_${timestamp}_stderr.log
   ```

## Source 04 Config (DO NOT LOSE)
These are source-specific overrides, not committed to main:
- **samples.json:** `{"samples": [{"Metadata_Source": "source_4"}]}`
- **sparcspy.yml:** nucleus diameter=25, cytosol diameter=69
- **profile/config.yaml:** 15 cores, vram_mb=80000, mem_mb=300000

## What Snakemake Will Do on Resume
- `--rerun-incomplete` will re-run any batches that were mid-extract when killed
- It will skip already-completed downloads/extractions (output files exist)
- Remaining work: ~19 downloads, ~19 extracts, then all 708 aggregate_broad_batch jobs, create_broad_structure, rechunk_broad

## Other Sources Still Pending
- source05: 81% done, source06: ~100% extraction, source07: 54%
- source09: 65%, source11: 74%, source13: ~100% extraction
- Each needs its own node; same resume pattern applies
