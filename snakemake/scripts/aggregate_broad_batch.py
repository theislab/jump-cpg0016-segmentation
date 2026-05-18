import json
import logging
import pathlib as pl

import h5py
import zarr

module_logger = logging.getLogger(__name__)


class StatusError(Exception):
    """Raised when segmentation/extraction status check fails"""

    pass


def _aggregate(
        extraction_path: str, segmentation_path: str, source: str, output: str
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

    output_base = pl.Path(output).joinpath("cellpose/objects")
    if not output_base.exists():
        msg = f"Output base directory ({output_base}) does not exist, aborting"
        module_logger.error(msg)
        raise FileNotFoundError(msg)

    with (
        h5py.File(segmentation_path, "r") as segmentation_file,
        h5py.File(extraction_path, mode="r") as extraction_file,
    ):
        cell_indices_by_image: dict[bytes, list[int]] = {}
        for row_idx, img_id in enumerate(extraction_file["single_cell_index_labelled"][:, 2]):
            cell_indices_by_image.setdefault(img_id, []).append(row_idx)

        zarr_groups: dict[str, zarr.Group] = {}
        mapping_written: set[pl.Path] = set()

        n_images = segmentation_file["labels"].shape[0]
        for image_idx in range(n_images):
            image_id = segmentation_file["labels"][image_idx, 1].decode()
            module_logger.info(f"Processing image {image_id}")
            _, batch, plate, _, _ = image_id.split("__")

            output_path = output_base.joinpath(batch, plate, f"{plate}.zarr")
            if not output_path.parent.exists():
                msg = (f"Output directory for zarr files does not exist, aborting. "
                       f"Batch {batch}, plate {plate}: ({output_path.parent})")
                module_logger.error(msg)
                raise FileNotFoundError(msg)

            if output_path.parent not in mapping_written:
                with output_path.parent.joinpath("channel_mapping.json").open("w") as mapping_file:
                    json.dump(obj=channel_mapping, fp=mapping_file)
                mapping_written.add(output_path.parent)

            cell_indices = cell_indices_by_image.get(image_id.encode(), [])

            store_key = str(output_path.resolve())
            if store_key not in zarr_groups:
                zarr_groups[store_key] = zarr.open(store_key, mode="a")
            out_zarr = zarr_groups[store_key]

            # overwrite=True makes resume after a TIMEOUT'd SLURM job idempotent: any
            # partial group left from a prior interrupted run gets wiped and rewritten
            # rather than raising ContainsGroupError. Each (plate, well, site) is
            # assigned to exactly one batch by design, so no concurrent writer can
            # collide on the same image_id.
            image_group = out_zarr.create_group(image_id, overwrite=True)
            image_group.create_array(name="label_image", data=segmentation_file["segmentation"][image_idx])
            image_group.create_array(name="single_cell_index", data=extraction_file["single_cell_index"][cell_indices])
            image_group.create_array(name="single_cell_data", data=extraction_file["single_cell_data"][cell_indices])

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
        extraction_path=snakemake.input.extraction_path,
        segmentation_path=snakemake.input.segmentation_path,
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
