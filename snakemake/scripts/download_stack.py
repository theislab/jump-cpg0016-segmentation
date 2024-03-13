import io
import logging
import multiprocessing
import pathlib as pl
import typing

import boto3
import botocore
import botocore.config as boto_conf
import numpy as np
import numpy.typing
import pandas as pd
from tifffile import tifffile

module_logger = logging.getLogger(__name__)

s3_client = boto3.client("s3", config=boto_conf.Config(signature_version=botocore.UNSIGNED))


class S3ERROR(Exception):
    """Raised when AWS S3 returns unexpected status codes"""

    pass


def get_bucket_key(s3_url: str) -> tuple[str, str]:
    """
    Gets bucket and key from a full s3 url

    Parameters
    ----------
    s3_url

    Returns
    -------
    Tuple of (bucket, key)

    """
    s3_url = pl.Path(s3_url)
    bucket = s3_url.parts[1]
    key = "/".join(s3_url.parts[2:])

    return bucket, key


def download_tif(s3_url: str) -> numpy.typing.NDArray[typing.Any]:
    """
    Downloads a tif to memory

    Parameters
    ----------
    s3_url
        URL of the tif file to download

    Returns
    -------
    np.ndarray of the image
    """
    module_logger.info(f"Downloading {s3_url}")
    bucket, key = get_bucket_key(s3_url)
    response = s3_client.get_object(Bucket=bucket, Key=key)

    if response["ResponseMetadata"]["HTTPStatusCode"] != 200:
        msg = f"Unexpected response from AWS s3: {response['ResponseMetadata']}"
        module_logger.error(msg)
        raise S3ERROR(msg)

    else:
        img = tifffile.imread(io.BytesIO(response["Body"].read()))
        return img


def download_stack(channel_urls: tuple[str]) -> numpy.typing.NDArray[typing.Any]:
    """
    Downloads and stacks

    Parameters
    ----------
    channel_urls
        OrderedDict of {channel name: s3 url}
    """
    images = [download_tif(url) for url in channel_urls]

    shape = None
    module_logger.info("Checking image shapes")
    for channel, image in enumerate(images):
        if shape is None and len(image.shape) == 2:
            shape = image.shape
        elif shape[0] == image.shape[0] and shape[1] == image.shape[1] and len(image.shape) == 2:
            continue
        else:
            msg = (
                f"Dimensions do not match up (shape: {image.shape}) for channel {channel},"
                f" s3 url {channel_urls[channel]}"
            )
            module_logger.critical(msg)
            raise AttributeError(msg)

    return np.stack(list(images))


def download_batch(batch_id: str, metadata: pd.DataFrame, output_path: str):
    """Downloads a single stack for an entire snakemake batch from the samples_metadata.parquet file"""
    metadata = metadata.loc[metadata["snakemake_batch"] == batch_id]
    stack_urls = [
        (row["s3_OrigDNA"], row["s3_OrigAGP"], row["s3_OrigER"], row["s3_OrigMito"], row["s3_OrigRNA"])
        for _, row in metadata.iterrows()
    ]

    # if not isinstance(n_processes, int):
    #     n_processes = multiprocessing.cpu_count()

    module_logger.info(f"Running download with 5 processes for batch {batch_id}")
    with multiprocessing.Pool(processes=5) as pool:
        stacks = pool.map(download_stack, stack_urls)

    shape = stacks[0].shape
    for stack_idx, stack in enumerate(stacks):
        for shape_idx, (expected, actual) in enumerate(zip(shape, stack.shape)):
            if expected != actual:
                module_logger.error(
                    f"Shape mismatch in stack {stack_idx}, shape {shape_idx}: Expected {expected}, got {actual}"
                )

    module_logger.info(f"Saving batch stack to {output_path}")
    np.save(file=output_path, arr=np.stack(stacks))


def main(snakemake):
    """
    Downloads and stacks all channels from a jump cellpainting sample

    Parameters
    ----------
    snakemake
        Magic object

    """
    # log_path = snakemake.log[0]
    # logging.critical(log_path)
    # logging.critical(type(log_path))
    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
        datefmt="%d.%m.%y %H:%M:%S",
    )

    module_logger.info("Starting download")
    download_batch(
        batch_id=snakemake.wildcards.batch,
        metadata=pd.read_parquet(snakemake.config["samples_meta"]),
        output_path=snakemake.output[0],
    )

    # urls = collections.OrderedDict(
    #     dna=snakemake.params.dna,
    #     rna=snakemake.params.rna,
    #     agp=snakemake.params.agp,
    #     er=snakemake.params.er,
    #     mito=snakemake.params.mito,
    # )
    # module_logger.info("URLs:")
    # for key, url in urls.items():
    #     module_logger.info(f"{key.upper()}: {url}")
    #
    # download_stack(
    #     channel_urls=urls,
    #     sample_path=snakemake.output[0],
    # )


if __name__ == "__main__":
    if "snakemake" in locals():
        main(locals()["snakemake"])
    else:
        logging.basicConfig(
            filename="../log/download/source_2__20210607_Batch_2__1053601879__K06__1.log",
            level=logging.DEBUG,
            format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
            datefmt="%d.%m.%y %H:%M:%S",
        )

        # test
        download_batch(
            batch_id="snakemake_batch_0",
            metadata=pd.read_parquet("../config/selected_metadata.parquet"),
            output_path="../results/images/source_2__1053600674.npy",
        )

        # sample_id = "source_8__J3__A1166127__A01__1"
        # sample_path = sample_id + ".tif"
        # channel_urls = collections.OrderedDict(
        #     dna="s3://cellpainting-gallery/cpg0016-jump/source_8/images/J3/images/A1166127/Images/HTS_A01_s1_w52C47F0FB-BD5E-465A-8A92-9E2ED9E72A12.tif",
        #     agp="s3://cellpainting-gallery/cpg0016-jump/source_8/images/J3/images/A1166127/Images/HTS_A01_s1_w29D9B15E0-0F42-4A59-A2F8-619BD017D1D1.tif",
        #     er="s3://cellpainting-gallery/cpg0016-jump/source_8/images/J3/images/A1166127/Images/HTS_A01_s1_w4194CE2FF-40C4-45F0-B9D7-AE81D0047EEA.tif",
        #     mito="s3://cellpainting-gallery/cpg0016-jump/source_8/images/J3/images/A1166127/Images/HTS_A01_s1_w1250BC2CE-02B0-48F6-9533-69E9CAC91BF5.tif",
        #     rna="s3://cellpainting-gallery/cpg0016-jump/source_8/images/J3/images/A1166127/Images/HTS_A01_s1_w30A163325-7B23-4E3C-BC17-E995BCF8818D.tif"
        # )
        # download_stack(
        #     # channel_urls=collections.OrderedDict(
        #     #     dna="s3://cellpainting-gallery/cpg0016-jump/source_8/images/J3/images/A1166127/Images/HTS_A01_s1_w52C47F0FB-BD5E-465A-8A92-9E2ED9E72A12.tif",
        #     #     agp="s3://cellpainting-gallery/cpg0016-jump/source_8/images/J3/images/A1166127/Images/HTS_A01_s1_w29D9B15E0-0F42-4A59-A2F8-619BD017D1D1.tif",
        #     #     er="s3://cellpainting-gallery/cpg0016-jump/source_8/images/J3/images/A1166127/Images/HTS_A01_s1_w4194CE2FF-40C4-45F0-B9D7-AE81D0047EEA.tif",
        #     #     mito="s3://cellpainting-gallery/cpg0016-jump/source_8/images/J3/images/A1166127/Images/HTS_A01_s1_w1250BC2CE-02B0-48F6-9533-69E9CAC91BF5.tif",
        #     #     rna="s3://cellpainting-gallery/cpg0016-jump/source_8/images/J3/images/A1166127/Images/HTS_A01_s1_w30A163325-7B23-4E3C-BC17-E995BCF8818D.tif",
        #     # ),
        #     channel_urls=collections.OrderedDict(
        #         dna="s3://cellpainting-gallery/cpg0016-jump/source_10/images/2021_08_17_U2OS_48_hr_run16/images/Dest210809-141456/Dest210809-141456_K06_T0001F001L01A01Z01C01.tif",
        #         agp="s3://cellpainting-gallery/cpg0016-jump/source_10/images/2021_08_17_U2OS_48_hr_run16/images/Dest210809-141456/Dest210809-141456_K06_T0001F001L01A03Z01C04.tif",
        #         er="s3://cellpainting-gallery/cpg0016-jump/source_10/images/2021_08_17_U2OS_48_hr_run16/images/Dest210809-141456/Dest210809-141456_K06_T0001F001L01A02Z01C02.tif",
        #         mito="s3://cellpainting-gallery/cpg0016-jump/source_10/images/2021_08_17_U2OS_48_hr_run16/images/Dest210809-141456/Dest210809-141456_K06_T0001F001L01A01Z01C05.tif",
        #         rna="s3://cellpainting-gallery/cpg0016-jump/source_10/images/2021_08_17_U2OS_48_hr_run16/images/Dest210809-141456/Dest210809-141456_K06_T0001F001L01A02Z01C03.tif",
        #     ),
        #     sample_path="source_10__2021_08_17_U2OS_48_hr_run16__Dest210809-141456__K06__1.tif",
        # )
