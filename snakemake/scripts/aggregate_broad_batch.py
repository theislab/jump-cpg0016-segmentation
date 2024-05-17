import json
import logging
import pathlib as pl

import h5py
import numpy as np
import zarr

module_logger = logging.getLogger(__name__)


class StatusError(Exception):
    """Raised when segmentation/extraction status check fails"""

    pass


def _aggregate(
    extraction_path: str, segmentation_path: str, source: str, output: str, debug: bool = False
) -> None:
    module_logger.info(f"Aggregating plates for {source}")
    channel_mapping = {
        0: "NucleusMask",
        1: "CellMask",
        2: "DNA",
        3: "AGP",
        4: "ER",
        5: "Mito",
        6: "RNA",
    }
    module_logger.info(f"Channel mapping: {channel_mapping}")

    # extraction_paths = [pl.Path(extraction_batch) for extraction_batch in extraction_batches]
    # segmentation_paths = [pl.Path(segmentation_batch) for segmentation_batch in segmentation_batches]

    output_base = pl.Path(output).joinpath("cellpose/objects")
    if not output_base.exists():
        msg = f"Output base directory ({output_base}) does not exist, aborting"
        module_logger.error(msg)
        raise FileNotFoundError(msg)
    # output_base.mkdir(parents=True, exist_ok=True)

    # for extraction_path, segmentation_path in zip(extraction_paths, segmentation_paths):
    module_logger.info(f"Processing snakemake batch {extraction_path.parent.name}")
    with (
        h5py.File(segmentation_path, "r") as segmentation_file,
        h5py.File(extraction_path, mode="r") as extraction_file,
    ):
        n_images = segmentation_file["labels"].shape[0]
        for image_idx in range(n_images):
            image_id = segmentation_file["labels"][image_idx, 1].decode()
            module_logger.info(f"Processing image {image_id}")
            source, batch, plate, well, spot = tuple(image_id.split("__"))

            output_path = output_base.joinpath(batch, plate, f"{plate}.zarr")
            if not output_path.parent.exists():
                msg = f"Output directory for zarr files in batch {batch}, plate {plate}: ({output_path.parent}) does not exist, aborting"
                module_logger.error(msg)
                raise FileNotFoundError(msg)

            with output_path.parent.joinpath("channel_mapping.json").open("w") as mapping_file:
                json.dump(obj=channel_mapping, fp=mapping_file)

            cell_indices = np.where(extraction_file["single_cell_index_labelled"][:, 2] == image_id.encode())[
                0
            ].tolist()

            with zarr.open(str(output_path.resolve()), mode="a") as out_zarr:
                image_group = out_zarr.create_group(image_id, overwrite=debug)
                image_group.array(name="label_image", data=segmentation_file["segmentation"][image_idx])
                image_group.array(name="single_cell_index", data=extraction_file["single_cell_index"][cell_indices])
                image_group.array(name="single_cell_data", data=extraction_file["single_cell_data"][cell_indices])

                module_logger.info(f"{image_id} wrote to {output_path.name}")


def main(snakemake):
    """
    Runs segmentation using SPARCSpy and CellPose

    Parameters
    ----------
    snakemake
        The "magic" snakemake object

    """
    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.DEBUG,
        format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
        datefmt="%d.%m.%y %H:%M:%S",
    )

    module_logger.info(f"Starting h5 aggregation for {snakemake.wildcards.source}")

    for key, value in snakemake.__dict__.items():
        module_logger.debug(f"Snakemake magic | {key}: {value}")

    _aggregate(
        extraction_path==snakemake.input.extraction_path,
        segmentation_path==snakemake.input.segmentation_path,
        source=snakemake.wildcards.source,
        output=snakemake.input.broad_dir,
    )

    module_logger.info("Aggregation complete")
    module_logger.info("Writing checkpoint file")

    with open(snakemake.output.checkpoint, "w") as checkpoint_file:
        checkpoint_file.write("Aggregation completed")


if __name__ == "__main__":
    if "snakemake" in locals():
        main(locals()["snakemake"])