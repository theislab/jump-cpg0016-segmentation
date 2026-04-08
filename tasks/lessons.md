# Lessons Learned

## System snakemake (v7) has broken pulp dependency
The system-installed snakemake 7.32.4 at `~/.local/bin/snakemake` fails with `AttributeError: module 'pulp' has no attribute 'list_solvers'` due to PuLP 2.9 removing that function. Fix: use pixi to install snakemake 8 (as source_03 does). Don't try to downgrade PuLP — the system has multiple Python versions (3.9 for snakemake, 3.12 for pip) making it messy.
**Why:** Lost 20 min debugging this. Use pixi for snakemake everywhere.

## pixi run needs cwd to contain pixi.toml
`pixi run` looks for `pixi.toml` in the current working directory. Use `--manifest-path` if running from elsewhere, or ensure the shell cd's first. In Claude Code, bash tool resets cwd between calls unless cd is chained in the same command.
**Why:** Multiple failed pipeline launch attempts due to wrong cwd.

## Snakemake directory locks persist across failed runs
If snakemake is killed or a concurrent instance tries to run, the `.snakemake` lock persists. Always `pixi run unlock` before re-launching. Check `pgrep -af snakemake` to ensure no stale processes hold the lock.
**Why:** Several launch attempts failed silently due to stale locks from previous attempts.

## zarr v3 API breaks: no context manager, `.array()` → `.create_array()`, no `codecs` kwarg
In zarr v3, `zarr.open()` returns a `Group` that does NOT support `with` (context manager). Use plain assignment instead. `group.array(name=..., data=...)` is removed — use `group.create_array(name=..., data=...)`. For sharded writes, use `create_array(chunks=..., shards=..., compressors=...)` — the old `codecs=[ShardingCodec(...).to_dict()]` parameter is gone.
**Why:** The conda env had zarr 3.1.5 installed. All three of these broke the aggregate and rechunk scripts simultaneously.

## pixi task `env.PATH` prepend shadows conda env Python
If `pixi.toml` sets `env = { PATH = "/path/to/miniforge3/bin:$PATH" }`, the base miniforge Python takes precedence over any conda env activated by snakemake via `source activate`. Fix: put miniforge3/bin at the END of PATH (`$PATH:/path/to/miniforge3/bin`). Conda is still findable but won't shadow the activated env's Python.
**Why:** Spent significant time debugging why `import dask` failed — snakemake was running base Python 3.12 (no dask) instead of the conda env's Python 3.14 (has dask). Debug output showing `sys.executable` was the key diagnostic.

## Concurrent zarr.open() in mode "a" races on store initialization
When using ProcessPoolExecutor, multiple workers calling `zarr.open(path, mode="a")` on a new store race to create zarr metadata. Fix: initialize the store once with `zarr.open(path, mode="w")` in the parent process before spawning workers. Workers then open with `mode="a"`.
**Why:** Got `ContainsGroupError` when 15 workers simultaneously tried to open a fresh store.

## Download .npy files disappearing doesn't mean download failed
The extract rule consumes the .npy but it's not marked `temp()` — however, if the pipeline is killed mid-run and restarted, snakemake may re-download. Missing .npy files with existing download logs means the pipeline was interrupted, not that downloads failed.
**Why:** Initially misdiagnosed as download script ending prematurely.
