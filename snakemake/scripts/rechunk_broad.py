import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import zarr
from numcodecs import blosc
from zarr.codecs import BloscCodec

module_logger = logging.getLogger(__name__)

SUBGROUPS = ["label_image", "single_cell_data", "single_cell_index"]
CELL_CHUNK_SIZE = 1
CELLS_PER_SHARD = 4096
COMPRESSOR = BloscCodec(cname="zstd", clevel=3, shuffle="shuffle")
BLOSC_THREADS = 1  # per-worker; each subprocess sets this independently


def _target_chunks_and_shard_shape(shape):
    if len(shape) == 0:
        return shape, None
    shard_chunk = (CELL_CHUNK_SIZE,) + shape[1:]
    target_chunks = (CELLS_PER_SHARD,) + shape[1:]
    return target_chunks, shard_chunk


def process_group(store_path, new_store_path, group_name):
    """Rechunk a single image group into per-cell sharded Zarr v3. Called in a subprocess."""
    blosc.set_nthreads(BLOSC_THREADS)

    src = zarr.open(store_path, mode="r")
    if group_name not in src:
        return

    src_group = src[group_name]
    dst = zarr.open(new_store_path, mode="a")

    for sub in SUBGROUPS:
        if sub not in src_group:
            continue
        src_ds = src_group[sub]
        data = np.asarray(src_ds)
        shard_shape, chunk_shape = _target_chunks_and_shard_shape(data.shape)
        if chunk_shape is None:
            module_logger.warning(f"Skipping {group_name}/{sub}: scalar array.")
            continue
        dst.create_array(
            name=f"{group_name}/{sub}",
            data=data,
            chunks=chunk_shape,
            shards=shard_shape,
            compressors=COMPRESSOR,
            overwrite=True,
        )


def rechunk_plate(store_path: Path, n_workers: int):
    new_store_path = Path(str(store_path).replace("/broad/", "/broad_compressed/"))
    os.makedirs(new_store_path.parent, exist_ok=True)

    src = zarr.open(str(store_path), mode="r")
    groups = list(src.group_keys())

    # Initialize the zarr store once before spawning workers
    zarr.open(str(new_store_path), mode="w")

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [
            executor.submit(process_group, str(store_path), str(new_store_path), group_name)
            for group_name in groups
        ]
        for future in as_completed(futures):
            future.result()  # raises if the subprocess failed


def main(snakemake):
    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
        datefmt="%d.%m.%y %H:%M:%S",
    )

    broad_dir = Path(snakemake.input.broad_dir)
    objects_dir = broad_dir / "cellpose/objects"
    n_workers = snakemake.threads

    module_logger.info(
        f"Rechunking plates under {objects_dir} with {n_workers} workers"
    )

    plate_paths = sorted(
        plate_dir / f"{plate_dir.name}.zarr"
        for batch_dir in sorted(objects_dir.iterdir()) if batch_dir.is_dir()
        for plate_dir in sorted(batch_dir.iterdir()) if plate_dir.is_dir()
        for _ in [None]
        if (plate_dir / f"{plate_dir.name}.zarr").exists()
    )

    module_logger.info(f"Found {len(plate_paths)} plates to rechunk")

    for plate_path in plate_paths:
        module_logger.info(f"Rechunking {plate_path.name}")
        rechunk_plate(plate_path, n_workers=n_workers)
        module_logger.info(f"Done: {plate_path.name}")

    with open(snakemake.output.checkpoint, "w") as f:
        f.write(f"Rechunked {len(plate_paths)} plates for {snakemake.wildcards.source}\n")

    module_logger.info("All plates rechunked successfully")


if __name__ == "__main__":
    if "snakemake" in locals():
        main(locals()["snakemake"])
