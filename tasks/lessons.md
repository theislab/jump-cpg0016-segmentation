# Lessons Learned

## conda activation fails in pixi task subprocesses

In a pixi task that sets `env = { PATH = "/home/icb/tim.treis/miniforge3/bin:$PATH" }`, snakemake's `conda:` + `script:` directive silently fails to activate the conda env. `source activate` runs but doesn't actually prepend the env to PATH — `which python` still returns miniforge3's Python 3.12 rather than the env's Python.

**Why:** The conda shell hook isn't initialized in a non-interactive bash subprocess when miniforge3 is already at the front of PATH.

**Fix:** Use `shell:` directive and invoke the conda env's Python directly by absolute path:
```smk
shell:
    """
    CONDA_PYTHON=/path/to/conda/envs/HASH/bin/python
    $CONDA_PYTHON scripts/myscript.py --arg {wildcards.x} >> {log} 2>&1
    """
```
This requires the script to have an argparse CLI (not just the magic `snakemake` object).

**Conda env hashes for this project:**
- Download env (boto3): `507f1bacc8e964624ea0c0fc2f042010_`
- Extract env (SPARCSpy): `79ff5f567d72279473b1703565acfd26_`
- Env prefix: `/ictstr01/groups/ml01/projects/2023_hackathon23_subcellular_spatial_niklas.schmacke/jump/envs/`

---

## pixi task `env.PATH` prepend shadows conda env Python

If `pixi.toml` sets `env = { PATH = "/path/to/miniforge3/bin:$PATH" }`, the base miniforge Python takes precedence over any conda env activated by snakemake via `source activate`. Fix: put miniforge3/bin at the END of PATH (`$PATH:/path/to/miniforge3/bin`). Conda is still findable but won't shadow the activated env's Python.
**Why:** Spent significant time debugging why `import dask` failed — snakemake was running base Python 3.12 (no dask) instead of the conda env's Python 3.14 (has dask). Debug output showing `sys.executable` was the key diagnostic.

---

## System snakemake (v7) has broken pulp dependency

The system-installed snakemake 7.32.4 at `~/.local/bin/snakemake` fails with `AttributeError: module 'pulp' has no attribute 'list_solvers'` due to PuLP 2.9 removing that function. Fix: use pixi to install snakemake 8 (as source_03 does). Don't try to downgrade PuLP — the system has multiple Python versions (3.9 for snakemake, 3.12 for pip) making it messy.
**Why:** Lost 20 min debugging this. Use pixi for snakemake everywhere.

---

## pixi run needs cwd to contain pixi.toml

`pixi run` looks for `pixi.toml` in the current working directory. Use `--manifest-path` if running from elsewhere, or ensure the shell cd's first. In Claude Code, bash tool resets cwd between calls unless cd is chained in the same command.
**Why:** Multiple failed pipeline launch attempts due to wrong cwd.

---

## Snakemake directory locks persist across failed runs

If snakemake is killed or a concurrent instance tries to run, the `.snakemake` lock persists. Always `pixi run unlock` before re-launching. Check `pgrep -af snakemake` to ensure no stale processes hold the lock.
**Why:** Several launch attempts failed silently due to stale locks from previous attempts.

---

## Lock files can be numbered — clear ALL of them

Snakemake creates `0.input.lock`, `0.output.lock`, `1.output.lock`, etc. Clearing only `0.*` misses numbered locks. Use:
```bash
find .snakemake/locks/ -name "*.lock" -delete
```
Not `rm -f .snakemake/locks/0.*.lock`.

---

## Snakemake gets stuck in shutdown on NFS after job failures

After all retries for a failed job are exhausted, snakemake prints "Shutting down, this might take some time" / "Exiting because a job execution failed" but the process never exits. It hangs indefinitely in Python multiprocessing cleanup blocked on NFS.

**Fix:** `kill -9 <PID>`, then clear stale locks before restarting:
```bash
rm -f .snakemake/locks/*.lock
```
`pixi run unlock` also works but rebuilds the full DAG (~30 min) — manual delete is faster.

---

## kill -9 on snakemake doesn't kill child processes

`kill -9` on the snakemake parent leaves all child bash/python processes as orphans. They continue running and consuming resources.

**Fix:** After killing snakemake, also kill orphaned children:
```bash
ps aux | grep "scripts/extract.py" | grep -v grep | awk '{print $2}' | xargs kill -9
```

---

## NFS write cache causes severe log lag

On this cluster, NFS write caching can cause log files and directory listings to appear frozen for 30–60+ minutes even when the pipeline is actively writing. `stat` on the log file will show a stale mtime. The process count (`ps aux | grep snakemake`) is the reliable liveness indicator.

**How to check real progress:** Count output files directly on disk rather than relying on log messages.

---

## Snakemake DAG build is slow on NFS with large job counts

For ~6000+ jobs, DAG construction takes 30–35 minutes on NFS. Plan for this when restarting the pipeline — the pipeline isn't hung, it's just building.

---

## GPU is NOT the bottleneck — NFS IO is

For this SPARCSpy cellpose pipeline, each batch has two phases:
- Cellpose segmentation (GPU): ~2 min per batch (fast)
- HDF5 cell extraction (CPU+NFS IO): ~16 min per batch (slow)

The GPU (A100 80GB) is idle ~85% of the time. The bottleneck is writing memory-mapped temp arrays and HDF5 files over NFS.

**What actually speeds things up:**
1. Local NVMe for `tmpdir` and extraction results — biggest win
2. More parallel jobs (raise `vram_mb`, `mem_mb`, `cores` in profile)
3. More CPU cores on the node
4. Bigger GPU does NOT help

**Working profile for 12 parallel jobs (confirmed stable):**
```yaml
cores: 12       # THE reliable concurrency cap — mem_mb doesn't work!
resources:
  vram_mb: 70000
  mem_mb: 240000
```
Note: `cores` is the only reliable way to limit concurrency. See "Snakemake resources doesn't reliably cap concurrency" lesson.

---

## Snakemake `resources` in profile doesn't reliably cap concurrency — use `cores`

Setting `resources: mem_mb: 240000` in the profile YAML with each job needing `mem_mb=20000` should cap at 12 jobs. In practice, snakemake 8 scheduled 33+ concurrent jobs, ignoring the constraint entirely.

**Fix:** Use `cores: 12` as a hard cap. Since extract jobs don't declare `threads`, each counts as 1 core. This reliably limits concurrency.

---

## zarr v3 API breaks: no context manager, `.array()` → `.create_array()`, no `codecs` kwarg

In zarr v3, `zarr.open()` returns a `Group` that does NOT support `with` (context manager). Use plain assignment instead. `group.array(name=..., data=...)` is removed — use `group.create_array(name=..., data=...)`. For sharded writes, use `create_array(chunks=..., shards=..., compressors=...)` — the old `codecs=[ShardingCodec(...).to_dict()]` parameter is gone.
**Why:** The conda env had zarr 3.1.5 installed. All three of these broke the aggregate and rechunk scripts simultaneously.

---

## Concurrent zarr.open() in mode "a" races on store initialization

When using ProcessPoolExecutor, multiple workers calling `zarr.open(path, mode="a")` on a new store race to create zarr metadata. Fix: initialize the store once with `zarr.open(path, mode="w")` in the parent process before spawning workers. Workers then open with `mode="a"`.
**Why:** Got `ContainsGroupError` when 15 workers simultaneously tried to open a fresh store.

---

## How to check per-source pipeline progress

Run `scripts/check_progress.sh` (see that file). It reports for each source:
- extractions done on disk vs total batches (from metadata parquet)
- images downloaded
- broad batches aggregated
- whether final aggregate exists

The key insight: **don't trust log file counts** (NFS lag). Count output files on disk instead.

---

## TiffFileError from S3 empty response is transient

`TiffFileError: not a TIFF file: header=b''` means S3 returned an empty body. This is a transient network issue. With `retries: 5` in the rule, the job gets 5 attempts. If all fail in one run, they get 5 fresh retries on the next restart.

---

## aggregate_broad.yml with ome-zarr installs zarr v2, not v3

The `aggregate_broad.yml` conda env spec originally depended on `ome-zarr`, which pulls in zarr 2.x. The updated `aggregate_broad_batch.py` uses `create_array()` which works in both v2 and v3, but `rechunk_broad.py` uses zarr v3-only features (`zarr.codecs.BloscCodec`, `shards=`, `group_keys()`). For consistency, change the dep to `zarr>=3` and drop `ome-zarr`.
**Why:** Discovered during source06 migration — the shared conda env (hash 425de6b...) had zarr 2.17.0. New env (hash d07e6420...) created with zarr v3.

---

## Per-source deployments diverge from repo main over time

Each `final_sourceNN/snakemake/` checkout drifts: local commits on main, stashed config changes, old branches. Migration pattern: stash → fetch → reset --hard origin/main → pop stash (resolve conflicts keeping main code + source-specific configs). The only files that should differ per-source are `config/samples.json`, `config/sparcspy.yml`, and `profile/config.yaml`.
**Why:** Found 6 different git states across 11 sources (two branches, two main versions, various local commits). Standardizing all at once prevents subtle differences in pipeline behavior.

---

## NFS cache coherency causes FileNotFoundError on newly created dirs

When `create_broad_structure` creates directories and `aggregate_broad_batch` jobs start immediately after (~7 min gap), some jobs get `FileNotFoundError` because NFS metadata cache hasn't propagated yet. Other parallel jobs accessing the same dirs succeed fine — it's a race.
**Why:** Source06 first aggregation run failed on batch_2389 with empty log (script didn't even start). All directories existed when checked manually afterward. Restart fixed it since dirs were fully materialized by then.
**Fix:** Just restart the pipeline. Snakemake resumes from completed checkpoints. Could also add a sleep between create_broad_structure and aggregate_broad_batch, but that requires rule changes.

---

## Git-committed config ≠ runtime config — verify extraction diameters by file dates, not git history

The sparcspy.yml in each deployment is a local override (not committed). Git history shows 21/53 diameters (852cc90, Mar 2024), but local files were updated to 23/57 in Apr 2026. To determine which diameters a given extraction used, check the **file modification timestamps** on extraction h5 files — not git log. If h5 files predate the local config change, they used old diameters and are invalid.

**Why:** Assumed all 2024 extractions used wrong diameters. This turned out to be correct, but source04 had been silently rerun in Apr 2026 with correct diameters. Checking only git history would have missed this. Always verify by timestamps.

---

## Snakemake replays full DAG even when only final step is missing

If intermediate outputs (extraction h5 files) have been deleted but their downstream checkpoints exist, snakemake still traces the full DAG and wants to rerun everything from scratch. `--force` on the target rule alone doesn't help — snakemake insists on recreating missing intermediates. If intermediates were intentionally deleted, the entire pipeline must rerun.

**Why:** Tried to run only rechunk_broad for source08 (broad aggregation existed, extraction deleted). Snakemake wanted 5988 jobs including all extractions. No shortcut — full rerun needed.

---

## Download .npy files disappearing doesn't mean download failed

The extract rule consumes the .npy but it's not marked `temp()` — however, if the pipeline is killed mid-run and restarted, snakemake may re-download. Missing .npy files with existing download logs means the pipeline was interrupted, not that downloads failed.
**Why:** Initially misdiagnosed as download script ending prematurely.
