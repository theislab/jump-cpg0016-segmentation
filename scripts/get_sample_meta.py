import logging

import pandas as pd

module_logger = logging.getLogger(__name__)


def download_df(url: str) -> pd.DataFrame:
    """
    Downloads a csv from url and reads it as df

    Parameters
    ----------
    url
        Of the csv file to download

    Returns
    -------
    pd.DataFrame
    """
    return pd.read_csv(url)


def _download_metadata(plate_link: str, well_link: str):
    module_logger.info(f"Downloading plate data from {plate_link}")
    plate = download_df(plate_link)
    module_logger.info(f"Downloading plate data from {well_link}")
    well = download_df(well_link)
    module_logger.info("Download completed")

    module_logger.info("Merging plate and well data")
    meta = pd.merge(plate, well, on=["Metadata_Source", "Metadata_Plate"])

    return meta


def _export_metadata(meta, meta_path):
    module_logger.info(f"Writing metadata to {meta_path}")
    meta.to_csv(meta_path)
    module_logger.info("Metadata download complete")


def _select_samples(meta: pd.DataFrame, **kwargs):
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
        if column in meta.columns:
            indices = meta.loc[meta[column] == value].index
            if len(indices) == 0:
                module_logger.info(f"Found no matching samples for column {column}, skipping!")
            else:
                if selected_indices is None:
                    selected_indices = set(indices)
                else:
                    selected_indices.intersection_update(set(indices))

        else:
            module_logger.warning(f"Column {column} not in meta dataframe, skipping!")

    selected_meta = meta.iloc[list(selected_indices)]
    return selected_meta


# def main():
#     # noinspection GrazieInspection
#     """Downloads metadata"""
#     logging.basicConfig(
#         level=logging.INFO,
#         format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
#         datefmt="%d.%m.%y %H:%M:%S",
#     )
#     meta = download_metadata(
#         plate_link="https://github.com/jump-cellpainting/datasets/raw/main/metadata/plate.csv.gz",
#         well_link="https://github.com/jump-cellpainting/datasets/raw/main/metadata/well.csv.gz",
#     )
#
#     try:
#         meta = select_samples(meta=meta, **whitelist)
#
#     except KeyError:
#         module_logger.info("Detected no whitelist, continuing with all metadata")
#
#     export_metadata(meta, meta_path=snakemake.output.meta)
#
#
# if __name__ == "__main__":
#     main()
