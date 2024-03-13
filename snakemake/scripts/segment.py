import logging
import pathlib as pl
import shutil

import numpy as np
from sparcscore.pipeline import extraction, project, workflows
from tifffile import tifffile

module_logger = logging.getLogger(__name__)


class NoCellsSegmentedError(Exception):
    """Raised when SPARCSpy segmentation does not find any cells"""

    pass


def _write_segmentation_status_file(path: pl.Path, success: bool = True):
    with path.open("wt") as f:
        if success:
            f.write("SEGMENTATION SUCCESSFUL")
        else:
            f.write("SEGMENTATION FAILED")


def _segment(input_path: str, project_dir: str, config_path: str, debug: bool = False):
    # module_logger.info(f"Segmentation for image {input_path}")
    # module_logger.info(f"Project directory: {project_dir}")

    sparcs_project = project.Project(
        location_path=project_dir,
        config_path=config_path,
        # segmentation_f=workflows.DAPISegmentationCellpose,
        segmentation_f=workflows.CytosolSegmentationCellpose,
        extraction_f=extraction.HDF5CellExtraction,
        debug=False,
        overwrite=True,
    )

    with tifffile.TiffFile(input_path) as tif:
        stack = np.array([page.asarray() for page in tif.pages])

    sparcs_project.load_input_from_array(stack)
    sparcs_project.segment()

    is_any_cell_segmented = np.any(sparcs_project.segmentation_f.maps["nucleus_segmentation"] > 0) and np.any(
        sparcs_project.segmentation_f.maps["cytosol_segmentation"] > 0
    )
    if is_any_cell_segmented:
        try:
            sparcs_project.extract()
        except:  # noqa: E722
            raise NoCellsSegmentedError  # noqa: B904
        finally:
            if not debug:
                try:
                    cache = pl.Path(sparcs_project.extraction_f.extraction_cache)
                    cache.rmdir()

                except AttributeError:
                    pass
    else:
        raise NoCellsSegmentedError


class NumbaLogFilter(logging.Filter):
    """Filters out numba debug logs"""

    def filter(self, record: logging.LogRecord) -> bool:
        """Filters out numba debug logs"""
        return record.levelno > logging.DEBUG or not record.name.startswith("numba")


def main(snakemake):
    """
    Runs segmentation using SPARCSpy and CellPose

    Parameters
    ----------
    snakemake
        The "magic" snakemake object

    """
    success_path = pl.Path(snakemake.output.success)
    extraction_target_path = pl.Path(snakemake.output.extraction)
    extraction_target_path.parent.mkdir(exist_ok=True, parents=True)

    try:
        _segment(
            input_path=snakemake.input.image,
            project_dir=snakemake.output.dir,
            config_path=snakemake.params.config,
            debug=snakemake.params.debug,
        )
        extraction_source_path = pl.Path(snakemake.output.dir).joinpath("extraction/data/single_cells.h5")
        shutil.move(extraction_source_path, extraction_target_path)

        _write_segmentation_status_file(success_path, success=True)

    except NoCellsSegmentedError:
        extraction_target_path.touch()
        _write_segmentation_status_file(success_path, success=False)


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

        sample = "source_10__2021_08_17_U2OS_48_hr_run16__Dest210809-141456__K06__3"
        # # Sample with no cells for testing error handling:
        # sample = "source_10__2021_08_03_U2OS_48_hr_run12__Dest210726-160150__L06__6"
        _segment(
            input_path=f"../results/images/{sample}.tif",
            project_dir=f"../results/segmentation/{sample}",
            config_path="../config/sparcspy.yml",
        )
