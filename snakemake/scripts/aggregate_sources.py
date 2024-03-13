import logging
import pathlib as pl

import h5py

module_logger = logging.getLogger(__name__)


class StatusError(Exception):
    """Raised when segmentation/extraction status check fails"""

    pass


def _aggregate(
    batches: list[str],
    source: str,
    output: str,
) -> None:
    module_logger.info(f"Aggregating h5 for {source}")
    # 'label_names'
    # 'single_cell_data'
    # 'single_cell_index'
    # 'single_cell_index_labelled'
    batch_paths = [pl.Path(batch) for batch in batches]
    output_path = pl.Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.is_file():
        module_logger.warning(f"Removing existing output file {output_path}")
        output_path.unlink()

    with h5py.File(output, "w-") as out_file:
        source_group = out_file.create_group(source)
        module_logger.info(f"Created source group: {source_group}")

        for batch_path in batch_paths:
            with h5py.File(batch_path, "r") as batch_file:
                sample_index_map = _get_sample_index_map(batch_file, source_group)
                _fill_out_file(batch_file, sample_index_map)

    # For debugging and inspection
    # import matplotlib.pyplot as plt
    # with h5py.File(output, "r") as out_file:
    #     h5_path = ("source_2/ZGRWVQNYTFGQLL-UHFFFAOYSA-N/"
    #                "source_2__20210607_Batch_2__1053601879__K06__1/single_cell_data")
    #     plt.imshow(out_file[h5_path][0, 4, ...])
    #     plt.show()


def _fill_out_file(batch_file: h5py.File, sample_index_map: dict[h5py.Group, list[int]]):
    dataset_names = ["single_cell_data", "single_cell_index"]
    for sample_group, cell_indices in sample_index_map.items():
        for dataset_name in dataset_names:
            shape = list(batch_file[dataset_name].shape)
            shape[0] = len(cell_indices)
            shape = tuple(shape)

            dataset = _get_dataset(dataset_name, group=sample_group, shape=shape)
            for idx, cell_idx in enumerate(cell_indices):
                dataset[idx] = batch_file[dataset_name][cell_idx]


def _get_sample_index_map(batch_file: h5py.File, source_group: h5py.Group) -> dict[h5py.Group, list[int]]:
    inchi_key_label_idx = _get_label_idx(batch_file, b"Metadata_InChIKey")
    id_label_idx = _get_label_idx(batch_file, b"id")
    n_cells = batch_file.get("single_cell_index_labelled").shape[0]
    sample_index_map = {}
    for cell_idx in range(n_cells):
        cell_inchi_key = _get_cell_label(batch_file, cell_idx, inchi_key_label_idx)
        module_logger.debug(f"Cell inchi key: {cell_inchi_key}")
        sample_id = _get_cell_label(batch_file, cell_idx, id_label_idx)
        try:
            inchi_group = _get_group(cell_inchi_key, source_group)
        except ValueError as err:
            module_logger.error(f"Could not create inchi group for cell {cell_idx} of sample {sample_id}")
            raise err
        sample_group = _get_group(sample_id, inchi_group)
        labels = [label.decode() for label in list(batch_file.get("label_names")[:])]
        for idx, label in enumerate(labels):
            if label not in ["cellid", "index"]:
                sample_group.attrs.create(name=label, data=_get_cell_label(batch_file, cell_idx, idx))

        try:
            sample_index_map[sample_group].append(cell_idx)
        except KeyError:
            sample_index_map[sample_group] = [cell_idx]

    return sample_index_map


def _get_label_idx(batch_file: h5py.File, label: bytes):
    return list(batch_file.get("label_names")[:]).index(label)


def _get_cell_label(file: h5py.File, cell_idx: int, label_idx: int):
    return file.get("single_cell_index_labelled")[cell_idx, label_idx].decode()


def _get_group(name: str, group: h5py.Group) -> h5py.Group:
    if name in group:
        return group[name]
    else:
        module_logger.info(f"Creating sub_group {name}")
        sub_group = group.create_group(name)
        return sub_group


def _get_dataset(name: str, group: h5py.Group, **kwargs) -> h5py.Dataset:
    if name in group:
        return group[name]
    else:
        dataset = group.create_dataset(name=name, **kwargs)
        module_logger.info(f"Created Dataset {dataset}")
        return dataset


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
        batches=snakemake.input.batches,
        source=snakemake.wildcards.source,
        output=snakemake.output.h5,
    )


if __name__ == "__main__":
    if "snakemake" in locals():
        main(locals()["snakemake"])
    else:
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
            datefmt="%d.%m.%y %H:%M:%S",
        )
        # input_paths = [
        #     "results/extraction/source_2/source_2__20210607_Batch_2__1053601879__K06__1/source_2__20210607_Batch_2__1053601879__K06__1.h5",
        #     "results/extraction/source_2/source_2__20210607_Batch_2__1053601879__K06__2/source_2__20210607_Batch_2__1053601879__K06__2.h5",
        #     "results/extraction/source_2/source_2__20210607_Batch_2__1053601879__K06__3/source_2__20210607_Batch_2__1053601879__K06__3.h5",
        #     "results/extraction/source_2/source_2__20210607_Batch_2__1053601879__K06__4/source_2__20210607_Batch_2__1053601879__K06__4.h5",
        #     "results/extraction/source_2/source_2__20210607_Batch_2__1053601879__K06__5/source_2__20210607_Batch_2__1053601879__K06__5.h5",
        #     "results/extraction/source_2/source_2__20210607_Batch_2__1053601879__K06__6/source_2__20210607_Batch_2__1053601879__K06__6.h5",
        # ]
        # status_paths = [
        #     "results/extraction/source_2/source_2__20210607_Batch_2__1053601879__K06__1/segmentation_status.txt",
        #     "results/extraction/source_2/source_2__20210607_Batch_2__1053601879__K06__2/segmentation_status.txt",
        #     "results/extraction/source_2/source_2__20210607_Batch_2__1053601879__K06__3/segmentation_status.txt",
        #     "results/extraction/source_2/source_2__20210607_Batch_2__1053601879__K06__4/segmentation_status.txt",
        #     "results/extraction/source_2/source_2__20210607_Batch_2__1053601879__K06__5/segmentation_status.txt",
        #     "results/extraction/source_2/source_2__20210607_Batch_2__1053601879__K06__6/segmentation_status.txt",
        # ]

        # input_paths = ["results/extraction/source_2/source_2__1053600674/extracted_single_cells.h5"]
        # input_paths = [
        #     "results/extraction/source_2/source_2__1053600674/extracted_single_cells.h5",
        #     "results/extraction/source_2/source_2__1053599503/extracted_single_cells.h5",
        #     "results/extraction/source_2/source_2__1053597936/extracted_single_cells.h5",
        #     "results/extraction/source_2/source_2__1086293911/extracted_single_cells.h5",
        #     "results/extraction/source_2/source_2__1086293492/extracted_single_cells.h5",
        #     "results/extraction/source_2/source_2__1086293133/extracted_single_cells.h5",
        #     "results/extraction/source_2/source_2__1086292884/extracted_single_cells.h5",
        #     "results/extraction/source_2/source_2__1086292389/extracted_single_cells.h5",
        #     "results/extraction/source_2/source_2__1086292037/extracted_single_cells.h5",
        #     "results/extraction/source_2/source_2__1086289686/extracted_single_cells.h5",
        # ]
        input_paths = [
            "results/extraction/source_5/source_5_snakemake_batch_0/extracted_single_cells.h5",
            "results/extraction/source_5/source_5_snakemake_batch_1/extracted_single_cells.h5",
            "results/extraction/source_5/source_5_snakemake_batch_2/extracted_single_cells.h5",
            "results/extraction/source_5/source_5_snakemake_batch_3/extracted_single_cells.h5",
            "results/extraction/source_5/source_5_snakemake_batch_4/extracted_single_cells.h5",
            "results/extraction/source_5/source_5_snakemake_batch_5/extracted_single_cells.h5",
            "results/extraction/source_5/source_5_snakemake_batch_6/extracted_single_cells.h5",
            "results/extraction/source_5/source_5_snakemake_batch_7/extracted_single_cells.h5",
            "results/extraction/source_5/source_5_snakemake_batch_8/extracted_single_cells.h5",
            "results/extraction/source_5/source_5_snakemake_batch_9/extracted_single_cells.h5",
            "results/extraction/source_5/source_5_snakemake_batch_10/extracted_single_cells.h5",
            "results/extraction/source_5/source_5_snakemake_batch_11/extracted_single_cells.h5",
            "results/extraction/source_5/source_5_snakemake_batch_12/extracted_single_cells.h5",
            "results/extraction/source_5/source_5_snakemake_batch_13/extracted_single_cells.h5",
        ]
        h5_path = "results/aggregated/source_5/source_5.h5"
        wildcard_source = "source_5"
        # inchi_key = "ZGRWVQNYTFGQLL-UHFFFAOYSA-N"
        _aggregate(
            batches=input_paths,
            source=wildcard_source,
            output=h5_path,
        )
