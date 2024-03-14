import importlib.util
import json
import logging
import multiprocessing as mp
import pathlib as pl

import boto3
import botocore
import botocore.config as boto_conf
import pandas as pd
import requests
from pandarallel import pandarallel
from tqdm import tqdm

engine_available = (
    importlib.util.find_spec("pyarrow") is not None or importlib.util.find_spec("fastparquet") is not None
)
if not engine_available:
    raise ImportError("Unable to find a usable engine; tried using: 'pyarrow', 'fastparquet'.")

module_logger = logging.getLogger(__name__)


class S3ERROR(Exception):
    """Raised when AWS S3 returns unexpected status codes"""

    pass


def parallel_s3_dl(paths: list[tuple[pl.Path, pl.Path]], n_processes: int = 10):
    """
    Downloads a list of s3 urls in parallel

    Parameters
    ----------
    paths
        tuple of s3 urls and local file paths to save the objects to
    n_processes
        number of processes to use for download
    """
    download_queue = mp.JoinableQueue()

    download_processes = [mp.Process(target=_download_worker, args=(download_queue,)) for _ in range(n_processes)]
    try:
        for process in download_processes:
            process.start()

        for path in paths:
            download_queue.put(path)

        download_queue.join()
        download_queue.close()

    finally:
        # noinspection PyUnresolvedReferences
        JumpDl.terminate_process_list(download_processes)


def _download_worker(download_queue: mp.JoinableQueue):
    client = boto3.client("s3", config=boto_conf.Config(signature_version=botocore.UNSIGNED))
    while True:
        s3_url, path = download_queue.get(block=True)

        bucket, key = JumpDl.get_bucket_key(s3_url=s3_url)
        response = client.get_object(Bucket=bucket, Key=key)

        if response["ResponseMetadata"]["HTTPStatusCode"] != 200:
            msg = f"Unexpected response from AWS s3: {response['ResponseMetadata']}"
            module_logger.error(msg)
            raise S3ERROR(msg)

        else:
            with path.open("wb") as f:
                f.write(response["Body"].read())

        download_queue.task_done()


class JumpConfig:
    """Base config class for jump downloader"""

    def __init__(self, jump_dir: pl.Path = None):
        self._setup_dirs(jump_dir)
        self._setup_paths()

        self._s3_client = None

        self._n_meta_processes = 7
        self._n_download_processes = 5

        self.profile_formatter = (
            "s3://cellpainting-gallery/cpg0016-jump/"
            "{Metadata_Source}/workspace/profiles/"
            "{Metadata_Batch}/{Metadata_Plate}/{Metadata_Plate}.parquet"
        )

        self.loaddata_formatter = (
            "s3://cellpainting-gallery/cpg0016-jump/"
            "{Metadata_Source}/workspace/load_data_csv/"
            "{Metadata_Batch}/{Metadata_Plate}/load_data_with_illum.parquet"
        )

        self._set_image_types()

    def _set_image_types(self):
        self._itypes = ["AGP", "DNA", "ER", "Mito", "RNA"]
        self.image_types = [f"Orig{itype}" for itype in self._itypes]

    def _setup_dirs(self, jump_dir):
        if isinstance(jump_dir, pl.Path):
            if jump_dir.is_dir():
                self.jump_dir = jump_dir
            else:
                module_logger.warning(f"Jump directory specified, but path is not a directory: {jump_dir}")
                module_logger.warning("Trying to create jump directory")
                jump_dir.mkdir(parents=True)
                self.jump_dir = jump_dir

        else:
            self.jump_dir = pl.Path.home().joinpath("Downloads/jump_cellpainting")
            module_logger.info(f"No jump directory specified, falling back to {self.jump_dir}")
            if not self.jump_dir.is_dir():
                module_logger.info("Trying to create jump directory")
                self.jump_dir.mkdir(parents=True)

        self.parquet_dir = self.jump_dir.joinpath("parquet")
        if not self.parquet_dir.is_dir():
            self.parquet_dir.mkdir()

        self.meta_dir = self.jump_dir.joinpath("metadata")
        if not self.meta_dir.is_dir():
            self.meta_dir.mkdir()

    def _setup_paths(self):
        self._plate_link = "https://github.com/jump-cellpainting/datasets/raw/main/metadata/plate.csv.gz"
        self._well_link = "https://github.com/jump-cellpainting/datasets/raw/main/metadata/well.csv.gz"
        self._compound_link = "https://github.com/jump-cellpainting/datasets/raw/main/metadata/compound.csv.gz"
        self._orf_link = "https://github.com/jump-cellpainting/datasets/raw/main/metadata/orf.csv.gz"

        self.plate_path = self.meta_dir.joinpath("plate.csv.gz")
        self.well_path = self.meta_dir.joinpath("well.csv.gz")
        self.compound_path = self.meta_dir.joinpath("compound.csv.gz")
        self.orf_path = self.meta_dir.joinpath("orf.csv.gz")

        self._meta_links = {
            self.plate_path: self._plate_link,
            self.well_path: self._well_link,
            self.compound_path: self._compound_link,
            self.orf_path: self._orf_link,
        }

        for path, link in self._meta_links.items():
            if not path.is_file():
                module_logger.info(f"Downloading {path}")
                file = requests.get(link)
                with path.open("wb") as f:
                    f.write(file.content)

    @staticmethod
    def get_bucket_key(s3_url: pl.Path) -> tuple[str, str]:
        """
        Gets bucket and key from a full s3 url

        Parameters
        ----------
        s3_url

        Returns
        -------
        Tuple of (bucket, key)

        """
        bucket = s3_url.parts[1]
        key = "/".join(s3_url.parts[2:])

        return bucket, key

    def s3_get(self, s3_url: pl.Path, path: pl.Path, client=None) -> None:
        """
        Downloads file from s3 url to a filesystem path

        Returns
        -------
        Nothing
        """
        bucket, key = self.get_bucket_key(s3_url)

        if client is None:
            if self._s3_client is None:
                client = boto3.client("s3", config=boto_conf.Config(signature_version=botocore.UNSIGNED))
            else:
                client = self._s3_client

        response = client.get_object(Bucket=bucket, Key=key)

        if response["ResponseMetadata"]["HTTPStatusCode"] != 200:
            msg = f"Unexpected response from AWS s3: {response['ResponseMetadata']}"
            module_logger.error(msg)
            raise S3ERROR(msg)

        else:
            with path.open("wb") as f:
                f.write(response["Body"].read())


class JumpMeta(JumpConfig):
    """Handles metadata for jump downloads"""

    def __init__(self, jump_dir: pl.Path = None):
        super().__init__(jump_dir)

        self._selected_indices = None
        self._selected_meta = None

        self._setup_dfs()

    def _setup_dfs(self):
        self.plate = pd.read_csv(self.plate_path)
        self.well = pd.read_csv(self.well_path)
        self.compound = pd.read_csv(self.compound_path)
        self.orf = pd.read_csv(self.orf_path)

        self.meta = pd.merge(self.plate, self.well, on=["Metadata_Source", "Metadata_Plate"])
        self.meta = pd.merge(self.meta, self.compound, on="Metadata_JCP2022", how="inner")
        self.meta.reset_index(drop=True)

    def whitelist_samples_from_json(self, path: pl.Path):
        # noinspection GrazieInspection
        """
        Selects samples from a json file

        JSON format:
        {
            "samples": [
                {
                    "column": "value",
                    "other column": "other value",
                    ...
                },
                {
                    # other sample
                }
            ]
        }

        Parameters
        ----------
        path: pl.Path to the json file

        Returns
        -------
        Nothing

        """
        module_logger.info(f"Using whitelist {path}")
        if not path.is_file() and path.suffix == ".json":
            msg = f"Sample file must be json, got: {path.absolute()}"
            module_logger.error(msg)
            raise ValueError(msg)
        else:
            with path.open("r") as f:
                samples = json.load(f)["samples"]

            for sample in samples:
                self.select_samples(**sample)

    def select_samples(self, **kwargs):
        """
        Iteratively select samples by matching column and values in meta table

        All kwargs are linked by boolean "and"

        Parameters
        ----------
        kwargs pairs of column names and values to select for
        """
        module_logger.debug(f"Selecting Sample: {kwargs}")
        selected_indices = None
        for column, value in kwargs.items():
            if column in self.meta.columns:
                indices = self.meta.loc[self.meta[column] == value].index
                if len(indices) == 0:
                    module_logger.info(f"Found no matching samples for column {column}, skipping!")
                else:
                    if selected_indices is None:
                        selected_indices = set(indices)
                    else:
                        selected_indices.intersection_update(set(indices))

            else:
                module_logger.warning(f"Column {column} not in meta dataframe, skipping!")

        if self._selected_indices is None:
            self._selected_indices = list(selected_indices)
        else:
            self._selected_indices.extend(selected_indices)

    def _get_selected_sample_data(self):
        if self._selected_indices is None:
            self._selected_meta = self.meta
            module_logger.info("Selecting all possible images")
        else:
            self._selected_indices = list(set(self._selected_indices))
            self._selected_meta = self.meta.iloc[self._selected_indices]

    def clear_samples(self):
        """
        Resets and clears the selected samples to None to start again with a new selection

        Returns
        -------
        Nothing
        """
        self._selected_indices = None
        self._selected_meta = None

    @staticmethod
    def id_from_row(row):
        """
        Retrieves source, batch, plate and well from a row of jump metadata

        Parameters
        ----------
        row pd.Series containing a row of jump metadata

        Returns
        -------
        source, batch, plate and well from a row of jump metadata
        """
        source = row["Metadata_Source"]
        batch = row["Metadata_Batch"]
        plate = row["Metadata_Plate"]
        well = row["Metadata_Well"]
        return source, batch, plate, well

    def path_from_id(self, source, batch, plate, well=None):
        """
        Builds a filesystem path relative to jump_dir for the specific plate/well

        Parameters
        ----------
        source: Metadata_Source
        batch: Metadata_Batch
        plate: Metadata_Plate
        well: Metadata_Well (optional)

        Returns
        -------
        pl.Path for the specific plate/well

        """
        if well is None:
            path = self.jump_dir.joinpath(source).joinpath(batch).joinpath(plate)
        else:
            path = self.jump_dir.joinpath(source).joinpath(batch).joinpath(plate).joinpath(well)
        path.mkdir(parents=True, exist_ok=True)
        return path


class JumpDl(JumpMeta):
    """Downloader class for JUMP cellpainting images and metadata"""

    def __init__(self, jump_dir: pl.Path = None, processes: int = None):
        super().__init__(jump_dir)
        if processes is not None:
            self._n_meta_processes = processes // 2
            self._n_download_processes = self._n_meta_processes

    def download_metadata(self):
        """
        Downloads selected samples

        Returns
        -------
        Nothing
        """
        self._get_dl_data()

    def download_samples(self):
        """
        Downloads selected samples

        Returns
        -------
        Nothing
        """

        self.download_metadata()

        meta_queue = mp.JoinableQueue()
        download_queue = mp.JoinableQueue()

        try:
            meta_processes = [
                mp.Process(target=_meta_worker, args=(meta_queue, download_queue))
                for _ in range(self._n_meta_processes)
            ]
            download_processes = [
                mp.Process(target=_jump_dl_worker, args=(download_queue,)) for _ in range(self._n_download_processes)
            ]

            for process in meta_processes:
                process.start()

            for process in download_processes:
                process.start()

            for _, row in tqdm(self._selected_meta.iterrows(), total=len(self._selected_meta)):
                meta_queue.put((self, row))

            meta_queue.join()
            meta_queue.close()

            download_queue.join()
            download_queue.close()

        finally:
            # noinspection PyUnboundLocalVariable
            self.terminate_process_list(meta_processes)
            # noinspection PyUnboundLocalVariable
            self.terminate_process_list(download_processes)

    @staticmethod
    def terminate_process_list(process_list: list[mp.Process]):
        """
        Terminates all processes in a list of processes

        Parameters
        ----------
        process_list
            List of mp.Process objects to terminate
        """
        for process in process_list:
            termination_counter = 0
            while process.is_alive():
                if termination_counter < 5:
                    process.terminate()
                else:
                    process.kill()

    def _get_dl_data(self):
        if self._selected_meta is None:
            self._get_selected_sample_data()

        load_data_index = ["Metadata_Source", "Metadata_Batch", "Metadata_Plate"]
        plates = self._selected_meta.set_index(load_data_index, drop=False).index.unique().to_list()

        download_paths = []
        parquet_paths = []
        for plate_id in tqdm(plates, total=len(plates)):
            parquet_path = self.parquet_dir.joinpath(f"{'__'.join(plate_id)}.parquet")
            parquet_paths.append(parquet_path)
            if not parquet_path.is_file():
                s3_path = self._get_parquet_url(*plate_id)
                download_paths.append((s3_path, parquet_path))

        parallel_s3_dl(download_paths, self._n_download_processes)

        load_data = pd.concat([pd.read_parquet(parquet_path) for parquet_path in parquet_paths])

        columns = ["Metadata_Source", "Metadata_Batch", "Metadata_Plate", "Metadata_Well"]
        self._selected_meta = self._selected_meta.merge(load_data, how="inner", on=columns)

        pandarallel.initialize(progress_bar=True)
        self._selected_meta = self._selected_meta.parallel_apply(self.add_id_col, axis=1)
        for channel in self.image_types:
            self._selected_meta = self._selected_meta.parallel_apply(self.add_s3_col, axis=1, channel=channel)

        snakemake_batch_column = "snakemake_batch"

        self._selected_meta[snakemake_batch_column] = generate_batch_ids(batch_size=250, metadata=self._selected_meta)

        columns = [
            "id",
            snakemake_batch_column,
            "Metadata_Source",
            "Metadata_PlateType",
            "Metadata_InChIKey",
            "Metadata_InChI",
        ]
        columns.extend([f"s3_{channel}" for channel in self.image_types])
        self._selected_meta = self._selected_meta[columns]

    @staticmethod
    def add_id_col(row: pd.Series):
        """
        Adds an id column

        Parameters
        ----------
        row
            _selected_meta row

        Returns
        -------
        row with "id" column, joining source, batch, plate, well and site using double underscores

        """
        row[
            "id"
        ] = f"{row['Metadata_Source']}__{row['Metadata_Batch']}__{row['Metadata_Plate']}__{row['Metadata_Well']}__{row['Metadata_Site']}"
        return row

    @staticmethod
    def add_s3_col(row: pd.Series, channel: str):
        """
        Adds joined path and filename for s3 urls

        Parameters
        ----------
        row
            _selected_meta row
        channel
            to create url for

        Returns
        -------
        _selected_meta row with the s3 url

        """
        row[f"s3_{channel}"] = row[f"PathName_{channel}"] + row[f"FileName_{channel}"]
        return row

    @staticmethod
    def get_image_url(row: pd.Series, image_type: str):
        """
        Builds the s3 url for an image based on a load_data row and the image type

        Parameters
        ----------
        row
            pd.Series row from load_data
        image_type
            str like "OrigAGP" etc.

        Returns
        -------
        tuple
            filename of the downloaded file and pl.Path with the s3 url
        """
        path = row[f"PathName_{image_type}"]
        file_name = row[f"FileName_{image_type}"]
        url = pl.Path(path).joinpath(file_name)
        return file_name, url

    @staticmethod
    def _get_parquet_url(source: str, batch: str, plate: str) -> pl.Path:
        url = (
            f"s3://cellpainting-gallery/cpg0016-jump/{source}/"
            f"workspace/load_data_csv/{batch}/{plate}/"
            f"load_data_with_illum.parquet"
        )
        return pl.Path(url)

    def export_meta(self, path: str = None):
        """
        Exports metadata of selected samples as csv

        Parameters
        ----------
        path
            Path of the metadata file

        Returns
        -------
        Nothing

        """
        try:
            if path is None:
                self._selected_meta.to_parquet(self.jump_dir.joinpath("selected_metadata.parquet"))
            else:
                self._selected_meta.to_parquet(path)
        except NameError:
            module_logger.error("No samples selected, please select samples first!")


def _meta_worker(meta_queue: mp.JoinableQueue, download_queue: mp.JoinableQueue):
    while True:
        jump_dl, row = meta_queue.get(block=True)
        assert isinstance(jump_dl, JumpDl)
        assert isinstance(row, pd.Series)

        jump_id = jump_dl.id_from_row(row)
        base_path = jump_dl.path_from_id(*jump_id)
        for image_type in jump_dl.image_types:
            file_name, image_url = jump_dl.get_image_url(row, image_type)
            image_site = "s" + row["Metadata_Site"]
            image_path = base_path.joinpath(f"{image_type}_{image_site}_{file_name}")
            if not image_path.is_file():
                download_queue.put((jump_dl, image_path, image_url))
        meta_queue.task_done()


def _jump_dl_worker(download_queue: mp.JoinableQueue):
    client = boto3.client("s3", config=boto_conf.Config(signature_version=botocore.UNSIGNED))
    while True:
        jump_dl, image_path, image_url = download_queue.get(block=True)
        assert isinstance(jump_dl, JumpDl)
        assert isinstance(image_path, pl.Path)
        assert isinstance(image_url, pl.Path)

        jump_dl.s3_get(image_url, image_path, client)
        download_queue.task_done()


def main(snakemake):
    # noinspection GrazieInspection
    """
    Main downloader in case this module is used in a snakemake workflow

    Parameters
    ----------
    snakemake
        Snakemake magic object

    Returns
    -------
    Nothing

    """
    logging.basicConfig(
        filename=snakemake.log.logfile,
        level=logging.INFO,
        format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
        datefmt="%d.%m.%y %H:%M:%S",
    )
    jdl = JumpDl(jump_dir=pl.Path(snakemake.output.jump_dir), processes=snakemake.threads)
    jdl.whitelist_samples_from_json(pl.Path(snakemake.input.samples))
    jdl.download_samples()
    jdl.export_meta(snakemake.output.metadata)
    module_logger.info("Download completed")


def generate_batch_ids(metadata: pd.DataFrame, batch_size: int = 250) -> list[str]:
    """Generates a list of batch ids based on Metadata_Source and a count (batch_size)"""
    source_index = metadata.columns.get_loc("Metadata_Source")
    batch_number = {}
    sample_counter = {}
    batch_text = "snakemake_batch"
    ids = []
    for idx in range(metadata.shape[0]):
        source = metadata.iat[idx, source_index]
        try:
            sample_counter[source] += 1
        except KeyError:
            sample_counter[source] = 1
            batch_number[source] = 0

        ids.append("_".join([source, batch_text, str(batch_number[source])]))
        if sample_counter[source] == batch_size:
            batch_number[source] += 1
            sample_counter[source] = 0

    return ids


if __name__ == "__main__":
    if "snakemake" in locals():
        main(locals()["snakemake"])