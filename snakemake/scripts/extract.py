import logging
import pathlib as pl
import shutil

import numpy as np
import pandas as pd
from sparcscore.pipeline import extraction, project, workflows


def get_metadata(metadata_path, batch: str):
    """Retrieves metadata for sparcspy batch processing from selected_metadata.parquet"""
    metadata = pd.read_parquet(metadata_path)
    batch_metadata = metadata.loc[metadata["snakemake_batch"] == batch]
    batch_metadata = batch_metadata.reset_index(drop=True)
    batch_metadata = batch_metadata.reset_index(drop=False)

    columns = [
        "index",
        "id",
        "snakemake_batch",
        "Metadata_Source",
        "Metadata_PlateType",
        "Metadata_InChIKey",
        "Metadata_InChI",
    ]

    return batch_metadata.get(columns)


# def _get_plate_image_stack(input_files: list[str], metadata: pd.DataFrame):
#     image_stacks = []
#     for image_path in input_files:
#         with tifffile.TiffFile(image_path) as tif:
#             image_stacks.append(np.array([page.asarray() for page in tif.pages]))
#
#     return np.stack(image_stacks)


def _find_single_tif_file(pattern: str, directory: pl.Path):
    file_matches = list(directory.glob(f"{pattern}.tif"))
    if len(file_matches) == 1:
        return file_matches[0]
    else:
        raise FileNotFoundError(f"Did not find unique file for id {pattern}")


def extraction_workflow(
    # input_files: list[str],
    batch_stack: str,
    project_dir: str,
    config_path: str,
    metadata_path: str,
    batch: str,
    debug: bool = False,
):
    """Runs batch segmentation and extraction using sparcspy"""
    sparcs_project = project.TimecourseProject(
        location_path=project_dir,
        config_path=config_path,
        # segmentation_f=workflows.DAPISegmentationCellpose,
        segmentation_f=workflows.Cytosol_Cellpose_BatchSegmentation,
        extraction_f=extraction.TimecourseHDF5CellExtraction,
        debug=debug,
        overwrite=False,
    )

    metadata = get_metadata(metadata_path, batch=batch)

    sparcs_project.load_input_from_array(np.load(batch_stack), label=metadata, overwrite=True)
    sparcs_project.segment()
    sparcs_project.extract()


def main(snakemake):
    """
    Runs segmentation using SPARCSpy and CellPose

    Parameters
    ----------
    snakemake
        The "magic" snakemake object

    """
    extraction_target_path = pl.Path(snakemake.output.extraction)
    extraction_target_path.parent.mkdir(exist_ok=True, parents=True)

    segmentation_target_path = pl.Path(snakemake.output.segmentation)
    segmentation_target_path.parent.mkdir(exist_ok=True, parents=True)

    # segmentation_target_path = extraction_target_path.parent.joinpath("segmentation.h5")

    extraction_workflow(
        batch_stack=snakemake.input.batch_stack,
        project_dir=snakemake.output.dir,
        config_path=snakemake.params.config,
        metadata_path=snakemake.config["samples_meta"],
        batch=snakemake.wildcards.batch,
        debug=snakemake.params.debug,
    )
    extraction_source_path = pl.Path(snakemake.output.dir).joinpath("extraction/data/single_cells.h5")
    segmentation_source_path = pl.Path(snakemake.output.dir).joinpath("segmentation/input_segmentation.h5")
    shutil.move(extraction_source_path, extraction_target_path)
    shutil.move(segmentation_source_path, segmentation_target_path)
    # shutil.move(segmentation_source_path, segmentation_target_path)

    #
    #     _write_segmentation_status_file(success_path, success=True)
    #
    # except NoCellsSegmentedError:
    #     extraction_target_path.touch()
    #     _write_segmentation_status_file(success_path, success=False)


class NumbaLogFilter(logging.Filter):
    """Filters out numba debug logs"""

    def filter(self, record: logging.LogRecord) -> bool:
        """Filters out numba debug logs"""
        return record.levelno > logging.DEBUG or not record.name.startswith("numba")


if __name__ == "__main__":
    if "snakemake" in locals():
        main(locals()["snakemake"])
    else:
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
            datefmt="%d.%m.%y %H:%M:%S",
        )
        for handler in logging.getLogger().handlers:
            handler.addFilter(NumbaLogFilter())

        extraction_workflow(
            batch_stack="../results/images/source_2__1053600674.npy",
            project_dir="../results/sparcspy/source_5/source_5_snakemake_batch_10/",
            config_path="../config/sparcspy.yml",
            metadata_path="../config/selected_metadata.parquet",
            batch="source_2__1053600674",
            debug=True,
        )

        # extraction_workflow(
        #     # input_files=[
        #     #     "../results/images/source_2__20210614_Batch_1__1053600674__G17__1.tif",
        #     #     "../results/images/source_2__20210614_Batch_1__1053600674__G17__2.tif",
        #     #     "../results/images/source_2__20210614_Batch_1__1053600674__G17__3.tif",
        #     #     "../results/images/source_2__20210614_Batch_1__1053600674__G17__4.tif",
        #     #     "../results/images/source_2__20210614_Batch_1__1053600674__G17__5.tif",
        #     #     "../results/images/source_2__20210614_Batch_1__1053600674__G17__6.tif",
        #     # ],
        #     input_files=[
        #         "../results/images/source_2__20210816_Batch_9__1086292389__G17__1.tif",
        #         "../results/images/source_2__20210816_Batch_9__1086292389__G17__2.tif",
        #         "../results/images/source_2__20210816_Batch_9__1086292389__G17__3.tif",
        #         "../results/images/source_2__20210816_Batch_9__1086292389__G17__4.tif",
        #         "../results/images/source_2__20210816_Batch_9__1086292389__G17__5.tif",
        #         "../results/images/source_2__20210816_Batch_9__1086292389__G17__6.tif",
        #     ],
        #     project_dir="../results/sparcspy/source_2/source_2__1053600674",
        #     config_path="../config/sparcspy.yml",
        #     # metadata_path="example_data/example_config_files/metadata_batch_cellpose.parquet",
        #     metadata_path="../config/selected_metadata.parquet",
        #     # batch="source_2__1053600674",
        #     batch="source_2__1086292389",
        #     debug=True,
        # )
