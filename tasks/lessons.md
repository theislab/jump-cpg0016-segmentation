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

## Download .npy files disappearing doesn't mean download failed
The extract rule consumes the .npy but it's not marked `temp()` — however, if the pipeline is killed mid-run and restarted, snakemake may re-download. Missing .npy files with existing download logs means the pipeline was interrupted, not that downloads failed.
**Why:** Initially misdiagnosed as download script ending prematurely.
