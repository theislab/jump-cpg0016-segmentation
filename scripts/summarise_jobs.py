"""Usage:
  summarise_jobs.py [--input_path=<input_path>] [--verbosity=<verbosity>] [--dry-run]

Arguments:
  -i --input_path=<input_path>    Path to the directory containing job files.
  -v --verbosity=<verbosity>      Verbosity level (0-3) [default: 1].

Options:
  -h --help                       Show this screen.
  --dry-run                       Perform a dry run without making any changes.
"""

from docopt import docopt
from pathlib import Path
from lamin_utils import logger
import sys

class PathValidationError(Exception):
    pass

def validate_paths(input_path: Path, verbosity: int, dry_run: bool):
    if not input_path:
        raise PathValidationError("Input path not provided.")
    
    if not input_path.exists():
        raise PathValidationError(f"The input path '{input_path}' does not exist.")
    
    if verbosity > 0:
        logger.info(f"Valid input path: {input_path}")

def summarise_jobs(input_path: Path, verbosity: int, dry_run: bool):
    job_files = list(input_path.glob("*.txt"))

    if verbosity > 0:
        logger.info(f"Found {len(job_files)} job file(s).")

    total_jobs = len(job_files)
    todo_jobs = 0
    done_jobs = 0
    error_jobs = 0

    for job_file in job_files:
        with open(job_file, "r") as f:
            status = f.readline().strip()
            if status == "todo":
                todo_jobs += 1
            elif status == "done":
                done_jobs += 1
            elif status == "error":
                error_jobs += 1

    if verbosity > 0:
        logger.info(f"Total jobs: {total_jobs}")
        logger.info(f"Todo jobs: {todo_jobs}")
        logger.info(f"Done jobs: {done_jobs}")
        logger.info(f"Error jobs: {error_jobs}")

    return total_jobs, todo_jobs, done_jobs, error_jobs

def main():
    arguments = docopt(__doc__, version="v0.1.0")

    input_path = arguments["--input_path"]
    dry_run = arguments["--dry-run"]
    verbosity = int(arguments["--verbosity"])

    input_path = Path(input_path).resolve() if input_path else None

    try:
        validate_paths(input_path, verbosity, dry_run)
    except PathValidationError as e:
        logger.error(e)
        sys.exit(1)

    if dry_run:
        logger.info("Dry run mode enabled. No changes will be made.")

    if verbosity > 0:
        logger.info(f"Dry run: {'enabled' if dry_run else 'disabled'}")
        logger.info(f"Verbosity level: {verbosity}")

    total_jobs, todo_jobs, done_jobs, error_jobs = summarise_jobs(input_path, verbosity, dry_run)

    if verbosity > 0:
        logger.success("Processing complete.")
        logger.info(f"Total jobs: {total_jobs}")
        logger.info(f"Todo jobs: {todo_jobs}")
        logger.info(f"Done jobs: {done_jobs}")
        logger.info(f"Error jobs: {error_jobs}")

if __name__ == "__main__":
    main()
