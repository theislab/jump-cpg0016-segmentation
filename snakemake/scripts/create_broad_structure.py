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


def _create_broad_structure(
    extraction_batches: list[str], segmentation_batches: list[str], source: str, output: str, debug: bool = False
) -> None:
    module_logger.info(f"Creating broad directory structure for {source}")

    extraction_paths = [pl.Path(extraction_batch) for extraction_batch in extraction_batches]
    segmentation_paths = [pl.Path(segmentation_batch) for segmentation_batch in segmentation_batches]

    output_base = pl.Path(output).joinpath("cellpose/objects")
    output_base.mkdir(parents=True, exist_ok=True)

    for extraction_path, segmentation_path in zip(extraction_paths, segmentation_paths):
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
                output_path.parent.mkdir(parents=True, exist_ok=True)

                # with output_path.parent.joinpath("channel_mapping.json").open("w") as mapping_file:
                #     json.dump(obj=channel_mapping, fp=mapping_file)

                # cell_indices = np.where(extraction_file["single_cell_index_labelled"][:, 2] == image_id.encode())[
                #     0
                # ].tolist()

                # with zarr.open(str(output_path.resolve()), mode="a") as out_zarr:
                #     image_group = out_zarr.create_group(image_id, overwrite=debug)
                #     image_group.array(name="label_image", data=segmentation_file["segmentation"][image_idx])
                #     image_group.array(name="single_cell_index", data=extraction_file["single_cell_index"][cell_indices])
                #     image_group.array(name="single_cell_data", data=extraction_file["single_cell_data"][cell_indices])

                #     module_logger.info(f"{image_id} wrote to {output_path.name}")


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

    _create_broad_structure(
        extraction_batches=snakemake.input.extraction_batches,
        segmentation_batches=snakemake.input.segmentation_batches,
        source=snakemake.wildcards.source,
        output=snakemake.output.broad_dir,
    )

    module_logger.info(f"Broad directory structure created for {snakemake.wildcards.source}")
    module_logger.info("Writing checkpoint file")

    with open(snakemake.output.checkpoint, "w") as checkpoint_file:
        checkpoint_file.write("Broad directory structure created for {snakemake.wildcards.source}")

if __name__ == "__main__":
    if "snakemake" in locals():
        main(locals()["snakemake"])
