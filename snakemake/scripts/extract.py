import logging
import pathlib as pl
import shutil

import h5py
import numpy as np
import pandas as pd
from sparcscore.pipeline import extraction, project, workflows


def _count_segmented_cells(project_dir: str) -> int:
    """Return number of cells that survived sparcspy's post-segmentation filter.

    Returns -1 if classes.csv is missing — callers should treat this as
    "unknown" and run the normal extract path so any real failure surfaces.
    """
    classes_csv = pl.Path(project_dir) / "segmentation" / "classes.csv"
    if not classes_csv.exists():
        return -1
    if classes_csv.stat().st_size == 0:
        return 0
    with classes_csv.open() as f:
        return sum(1 for _ in f)


def _write_empty_extraction_h5(target_path: pl.Path) -> None:
    """Write a schema-valid extraction h5 with zero rows.

    Why: when sparcspy's nucleus/cytosol pairing rejects every cell in a batch,
    sparcs_project.extract() crashes inside h5py with
    `IndexError: Index (0) out of range for empty dimension`, killing the whole
    snakemake DAG. Mirroring the schema with 0-row datasets lets
    aggregate_broad_batch.py iterate to zero per-image cells and move on.
    """
    label_names = [
        b"index", b"cellid", b"id", b"snakemake_batch",
        b"Metadata_Source", b"Metadata_PlateType",
        b"Metadata_InChIKey", b"Metadata_InChI",
    ]
    target_path.parent.mkdir(exist_ok=True, parents=True)
    vlen_bytes = h5py.special_dtype(vlen=bytes)
    with h5py.File(target_path, "w") as f:
        f.create_dataset("label_names", data=np.array(label_names, dtype=object), dtype=vlen_bytes)
        # No chunks/compression on the empty dataset: chunk shape can't exceed
        # data shape in any dimension, and lzf+chunks gain nothing on 0 rows.
        # aggregate_broad_batch reads by row index, so layout doesn't matter.
        f.create_dataset("single_cell_data", shape=(0, 7, 150, 150), dtype=np.float16)
        f.create_dataset("single_cell_index", shape=(0, 2), dtype=np.uint64)
        f.create_dataset("single_cell_index_labelled", shape=(0, 8), dtype=vlen_bytes)


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

    n_cells = _count_segmented_cells(project_dir)
    if n_cells == 0:
        logging.warning(
            "Batch %s produced 0 segmented cells; writing empty extraction h5 "
            "instead of running sparcs_project.extract() (which would crash on "
            "an empty single_cell_data shape).",
            batch,
        )
        empty_target = pl.Path(project_dir) / "extraction" / "data" / "single_cells.h5"
        _write_empty_extraction_h5(empty_target)
    else:
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
    if "snakemake" in dir():
        main(snakemake)  # noqa: F821 — injected by snakemake preamble
    else:
        import argparse

        parser = argparse.ArgumentParser(description="Run SPARCSpy extraction workflow")
        parser.add_argument("--batch-stack", required=True)
        parser.add_argument("--project-dir", required=True)
        parser.add_argument("--config", required=True)
        parser.add_argument("--metadata", required=True)
        parser.add_argument("--batch", required=True)
        parser.add_argument("--extraction", required=True)
        parser.add_argument("--segmentation", required=True)
        parser.add_argument("--debug", action="store_true", default=False)
        args = parser.parse_args()

        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
            datefmt="%d.%m.%y %H:%M:%S",
        )
        for handler in logging.getLogger().handlers:
            handler.addFilter(NumbaLogFilter())

        extraction_target_path = pl.Path(args.extraction)
        extraction_target_path.parent.mkdir(exist_ok=True, parents=True)
        segmentation_target_path = pl.Path(args.segmentation)
        segmentation_target_path.parent.mkdir(exist_ok=True, parents=True)

        extraction_workflow(
            batch_stack=args.batch_stack,
            project_dir=args.project_dir,
            config_path=args.config,
            metadata_path=args.metadata,
            batch=args.batch,
            debug=args.debug,
        )
        extraction_source_path = pl.Path(args.project_dir).joinpath("extraction/data/single_cells.h5")
        segmentation_source_path = pl.Path(args.project_dir).joinpath("segmentation/input_segmentation.h5")
        shutil.move(extraction_source_path, extraction_target_path)
        shutil.move(segmentation_source_path, segmentation_target_path)

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
