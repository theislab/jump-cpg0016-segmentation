"""Usage:
  transfer_plates.py [--input_path=<input_path>] [--verbosity=<verbosity>] [--dry-run]

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
import os
import json
import subprocess
import configparser
import boto3
import time
from datetime import datetime
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import re

S3_TARGET_DIR = "TIM2"


class PathValidationError(Exception):
    pass

def generate_target_prefix(local_path: str):
    """
    Generate the target prefix for S3.

    Args:
    local_path (str): The local path of the file to be uploaded.

    Returns:
    str: The target prefix for the file on S3.
    """

    # 1) extract the path in a source-num robust way
    # Regex explanation:
    #  - cpg\d+-jump  matches 'cpg' followed by one or more digits, then '-jump'
    #  - /source_\d+  matches '/source_' followed by one or more digits
    #  - (.*)$        captures the rest of the path until the end
    pattern = re.compile(r'(cpg0016-jump\/source_\d+.*)$')

    match = pattern.search(local_path)
    return f"{S3_TARGET_DIR}/{match[1]}" if match else None

def validate_paths(input_path: Path, verbosity: int):
    if not input_path:
        raise PathValidationError("Input path not provided.")

    if not input_path.exists():
        raise PathValidationError(f"The input path '{input_path}' does not exist.")

    if verbosity > 0:
        logger.info(f"Valid input path: {input_path}")


def log_into_job_file(job_file: Path, message: str):
    """
    Append a log message to the job file.

    Args:
    job_file (Path): The path to the job file.
    message (str): The log message to append to the job file.
    """
    try:
        timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        with open(job_file, "a") as f:
            f.write(f"[{timestamp}] {message}\n")
    except Exception as e:
        logger.error(f"Failed to append message to '{job_file.name}': {e}")


def update_job_status(job_file: Path, new_status: str, verbosity: int = 1):
    """
    Update the status of the job in the job file.

    Args:
    job_file (Path): The path to the job file.
    new_status (str): The new status to be written to the job file.
    """
    try:
        with open(job_file, "r+") as f:
            lines = f.readlines()
            if lines:
                lines[0] = f"{new_status}\n"  # Update the status line
            f.seek(0)
            f.writelines(lines)
            f.truncate()
        if verbosity > 1:
            logger.info(f"Updated status of '{job_file.name}' to '{new_status}'.")
    except Exception as e:
        if verbosity > 1:
            msg = f"Failed to update status of '{job_file.name}' to '{new_status}': {e}"
            logger.error(msg)
            log_into_job_file(job_file, msg)


def read_job_file(job_file: Path):
    with open(job_file, "r") as f:
        status = f.readline().strip()
        path = f.readline().strip()
        end_marker = f.readline().strip()
        assert end_marker == "---", "Invalid job file format."
        return status, path


def find_plate_to_transfer(input_folder: Path, verbosity: int):

    if not input_folder.exists() or not input_folder.is_dir():
        raise ValueError(
            f"The input folder '{input_folder}' does not exist or is not a directory."
        )

    job_files = sorted(input_folder.glob("*.txt"))

    if verbosity > 0:
        logger.info(f"Found {len(job_files)} job file(s) in '{input_folder}'.")
        logger.info("Scanning for 'todo' status...")

    # scan files until we find one with the status "todo"
    for job_file in job_files:
        with open(job_file, "r") as f:
            status = f.readline().strip()
            if status == "todo":
                if verbosity > 0:
                    logger.info(f"Found 'todo' job: {job_file.name}")
                return job_file

    if verbosity > 0:
        logger.warning("No 'todo' jobs found.")

    return None


def load_aws_credentials(profile_name="default"):
    """Load AWS credentials from the ~/.aws/credentials file."""
    config = configparser.ConfigParser()
    credentials_path = os.path.join(os.path.expanduser("~"), ".aws", "credentials")
    try:
        config.read(credentials_path)
        access_key = config.get(profile_name, "aws_access_key_id")
        secret_key = config.get(profile_name, "aws_secret_access_key")
        return access_key, secret_key
    except configparser.Error as e:
        logger.error(f"Error reading AWS credentials: {e}")
        raise


def get_temporary_credentials(custom_env):
    """Retrieve temporary AWS credentials using an AWS CLI command."""
    command = [
        "aws",
        "s3control",
        "get-data-access",
        "--account-id",
        "309624411020",
        "--target",
        "s3://staging-cellpainting-gallery/*",
        "--permission",
        "READWRITE",
        "--privilege",
        "Default",
        "--duration-seconds",
        "14400",  # Setting to 4 h
        "--region",
        custom_env["AWS_REGION"],
    ]
    try:
        result = subprocess.run(
            command, stdout=subprocess.PIPE, text=True, env=custom_env
        )
        result.check_returncode()  # Will raise CalledProcessError if the command failed
        data = json.loads(result.stdout)
        return data.get("Credentials")
    except subprocess.CalledProcessError as e:
        logger.error(f"CLI command failed: {e}")
        raise


def initialize_s3_client(credentials, region):
    """Initialize and return a Boto3 S3 client."""
    return boto3.client(
        "s3",
        aws_access_key_id=credentials["AccessKeyId"],
        aws_secret_access_key=credentials["SecretAccessKey"],
        aws_session_token=credentials["SessionToken"],
        region_name=region,
    )


def upload_file(s3_client, bucket_name, local_path, s3_path, job_file):
    """Upload a single file to an S3 bucket with error handling."""
    try:
        s3_client.upload_file(local_path, bucket_name, s3_path)
    except Exception as e:
        msg = f"Failed to upload {local_path} to {s3_path}: {e}"
        logger.error(msg)
        log_into_job_file(job_file, msg)
        return False
    return True


def upload_files(s3_client, bucket_name, local_directory, s3_directory, job_file):
    """Upload files to an S3 bucket using ThreadPoolExecutor for parallel processing."""
    files_to_upload = []
    for root, dirs, files in os.walk(local_directory):
        for file in files:
            local_path = os.path.join(root, file)
            s3_path = generate_target_prefix(local_path)
            files_to_upload.append((local_path, s3_path))

    try:
        with ThreadPoolExecutor() as executor:
            futures = [
                executor.submit(
                    upload_file, s3_client, bucket_name, local_path, s3_path, job_file
                )
                for local_path, s3_path in files_to_upload
            ]
            with tqdm(total=len(futures), desc="Uploading Files") as tqdm_instance:
                _upload_in_parallel(futures, tqdm_instance, job_file)
    except KeyboardInterrupt:
        error_message = "Transfer interrupted by user."
        log_into_job_file(job_file, error_message)
        raise

def _upload_in_parallel(futures, tqdm_instance, job_file):
    try:
        completed = 0
        for future in as_completed(futures):
            if future.result():
                tqdm_instance.update(1)
                completed += 1
                if completed % 10000 == 0:
                    progress_message = (
                        f"{completed}/{len(futures)} files uploaded."
                    )
                    log_into_job_file(job_file, progress_message)
        # Final log after last file
        progress_message = f"{completed}/{len(futures)} files uploaded."
        log_into_job_file(job_file, progress_message)
        msg = f"Successfully uploaded {completed} files."
        logger.success(msg)
        log_into_job_file(job_file, msg)
        update_job_status(job_file, "done")
    except Exception as e:
        msg = f"An error occurred: {e}"
        logger.error(msg)
        log_into_job_file(job_file, msg)
        raise


def main():
    arguments = docopt(__doc__, version="v0.1.0")

    input_path = arguments["--input_path"]
    dry_run = arguments["--dry-run"]
    verbosity = int(arguments["--verbosity"])

    input_path = Path(input_path).resolve() if input_path else None

    try:
        validate_paths(input_path, verbosity)
    except PathValidationError as e:
        logger.error(e)
        sys.exit(1)

    if dry_run:
        logger.info("Dry run mode enabled. No changes will be made.")

    if verbosity > 0:
        logger.info(f"Dry run: {'enabled' if dry_run else 'disabled'}")
        logger.info(f"Verbosity level: {verbosity}")

    while True:
        # Find the first job with status "todo"
        plate_to_transfer = find_plate_to_transfer(input_path, verbosity)
        if not plate_to_transfer:
            logger.info("No 'todo' jobs found. Exiting.")
            sys.exit(0)

        _, plate_path = read_job_file(plate_to_transfer)

        try:

            access_key_id, secret_access_key = load_aws_credentials("cpg_staging")
            custom_env = os.environ.copy()
            custom_env.update(
                {
                    "AWS_ACCESS_KEY_ID": access_key_id,
                    "AWS_SECRET_ACCESS_KEY": secret_access_key,
                    "AWS_SESSION_TOKEN": "",  # reset session token
                    "AWS_REGION": "us-east-1",
                }
            )

            if credentials := get_temporary_credentials(custom_env):
                if s3_client := initialize_s3_client(
                    credentials, custom_env["AWS_REGION"]
                ):
                    update_job_status(plate_to_transfer, "ongoing")

                    if not dry_run:
                        logger.info(f"Uploading files from '{plate_path}' to S3...")
                        if upload_status := upload_files(
                            s3_client,
                            "staging-cellpainting-gallery",
                            plate_path,
                            S3_TARGET_DIR,
                            plate_to_transfer,
                        ):
                            update_job_status(plate_to_transfer, "success")
                        else:
                            update_job_status(plate_to_transfer, "error")
                    else:
                        logger.info("Dry run: No files uploaded.")
                        update_job_status(plate_to_transfer, "todo")

            else:
                logger.error("No credentials retrieved.")
                update_job_status(plate_to_transfer, "error")

        except KeyboardInterrupt:
            msg = "Keyboard interrupt received. Exiting." 
            logger.error(msg)
            log_into_job_file(plate_to_transfer, msg)
            if plate_to_transfer:
                update_job_status(plate_to_transfer, "error")
            sys.exit(1)

        except Exception as e:
            msg = f"An error occurred: {e}"
            logger.error(msg)
            log_into_job_file(plate_to_transfer, msg)
            if plate_to_transfer:
                update_job_status(plate_to_transfer, "error")
            time.sleep(3)  # wait before retrying


if __name__ == "__main__":
    main()
