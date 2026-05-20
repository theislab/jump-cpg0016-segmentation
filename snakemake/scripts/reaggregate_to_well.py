"""Reaggregate Option A per-FOV zarr layout to Option B per-well layout.

Schema matches the Option B spec in the GitHub issue (well-level grouping for Broad
delivery; full source__batch__plate__well group names; src_row/cellid carried
verbatim from Option A; FOV identity preserved via a separate fov_index array).

Input  (Option A, per-FOV groups under one {plate}.zarr):
  {plate}.zarr/
    {source}__{batch}__{plate}__{well}__{fov}/
      label_image           uint32 (2, H, W)
      single_cell_index     uint64 (N, 2)         cols = [src_row, cellid]
      single_cell_data      f16    (N, 7, 150, 150)

Output (Option B, per-well groups under {plate}.zarr in sibling tree):
  {plate}.zarr/
    attrs: source, batch, plate, fov_height, fov_width, n_wells, created_at
    {source}__{batch}__{plate}__{well}/         # fully qualified group name
      attrs: well, n_fovs, n_cells
      label_image           uint32 (n_fovs, 2, H, W)  sharded per well
      single_cell_data      f16    (n_cells, 7, 150, 150)
      single_cell_index     uint64 (n_cells, 2)  cols = [src_row, cellid]  (verbatim)
      fov_index             uint32 (n_cells,)    1-indexed FOV number per cell
  channel_mapping.json (sibling, copied verbatim)

Failure model:
  - Refuses to start if destination {plate}.zarr already exists.
  - Writes to {plate}.zarr.tmp/; renames at the very end (atomic on Lustre).
  - Any leftover .tmp from a prior crash is wiped before retry.
  - Validates a sample of cells via (fov_index, cellid) -> mask lookup before commit.
"""

import argparse
import logging
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import zarr
from zarr.codecs import BloscCodec

module_logger = logging.getLogger(__name__)

SHARD_CELLS = 4096
COMPRESSOR = BloscCodec(cname="zstd", clevel=3, shuffle="shuffle")


def parse_group_name(name: str) -> tuple[str, str, str, str, int]:
    parts = name.split("__")
    if len(parts) != 5:
        raise ValueError(f"Unexpected group name (expect 5 '__'-separated parts): {name}")
    source, batch, plate, well, fov = parts
    return source, batch, plate, well, int(fov)


def reaggregate_plate(
    src_plate_zarr: Path,
    dst_plate_zarr: Path,
    wells_filter: list[str] | None = None,
    validate_sample: int = 5,
) -> None:
    module_logger.info(f"src: {src_plate_zarr}")
    module_logger.info(f"dst: {dst_plate_zarr}")

    if dst_plate_zarr.exists():
        raise RuntimeError(
            f"Destination already exists, refusing to overwrite: {dst_plate_zarr}. "
            f"Delete it manually if you want to re-run."
        )

    tmp_path = dst_plate_zarr.parent / f"{dst_plate_zarr.name}.tmp"
    if tmp_path.exists():
        module_logger.info(f"Removing leftover tmp from prior crash: {tmp_path}")
        shutil.rmtree(tmp_path)

    src = zarr.open(str(src_plate_zarr), mode="r")

    well_to_fovs: dict[str, list[tuple[int, str]]] = {}
    plate_id: tuple[str, str, str] | None = None
    for group_name in src.group_keys():
        source, batch, plate, well, fov_int = parse_group_name(group_name)
        if plate_id is None:
            plate_id = (source, batch, plate)
        elif (source, batch, plate) != plate_id:
            raise RuntimeError(
                f"Mixed plate IDs in {src_plate_zarr}: expected {plate_id}, "
                f"got group {group_name}"
            )
        if wells_filter and well not in wells_filter:
            continue
        well_to_fovs.setdefault(well, []).append((fov_int, group_name))

    if not well_to_fovs:
        raise RuntimeError(f"No FOV groups found in {src_plate_zarr} (after filter)")
    for w in well_to_fovs:
        well_to_fovs[w].sort()

    plate_source, plate_batch, plate_name = plate_id  # type: ignore[misc]
    module_logger.info(
        f"plate ({plate_source}, {plate_batch}, {plate_name}): "
        f"{len(well_to_fovs)} wells, {sum(len(v) for v in well_to_fovs.values())} FOVs"
    )

    first_well = next(iter(well_to_fovs))
    _, first_group = well_to_fovs[first_well][0]
    sample_li = src[first_group]["label_image"]
    if sample_li.ndim != 3 or sample_li.shape[0] != 2:
        raise RuntimeError(f"Unexpected label_image shape: {sample_li.shape}")
    fov_height = int(sample_li.shape[-2])
    fov_width = int(sample_li.shape[-1])

    dst_plate_zarr.parent.mkdir(parents=True, exist_ok=True)
    dst = zarr.open(str(tmp_path), mode="w")
    dst.attrs["source"] = plate_source
    dst.attrs["batch"] = plate_batch
    dst.attrs["plate"] = plate_name
    dst.attrs["fov_height"] = fov_height
    dst.attrs["fov_width"] = fov_width
    dst.attrs["n_wells"] = len(well_to_fovs)
    dst.attrs["created_at"] = datetime.now(timezone.utc).isoformat()

    for well, fovs in well_to_fovs.items():
        module_logger.info(f"  well {well}: {len(fovs)} FOVs")
        _write_well(
            src, dst,
            plate_source, plate_batch, plate_name,
            well, fovs, fov_height, fov_width,
        )

    src_cm = src_plate_zarr.parent / "channel_mapping.json"
    if src_cm.exists():
        dst_cm = dst_plate_zarr.parent / "channel_mapping.json"
        shutil.copyfile(src_cm, dst_cm)
        module_logger.info(f"copied {src_cm.name} -> {dst_cm}")
    else:
        module_logger.warning(f"missing channel_mapping.json: {src_cm}")

    _validate_round_trip(dst, validate_sample)

    tmp_path.rename(dst_plate_zarr)
    module_logger.info(f"committed: {dst_plate_zarr}")


def _write_well(
    src: zarr.Group,
    dst: zarr.Group,
    plate_source: str,
    plate_batch: str,
    plate_name: str,
    well: str,
    fovs: list[tuple[int, str]],
    fov_height: int,
    fov_width: int,
) -> None:
    n_fovs = len(fovs)
    fov_group_names = [gn for _, gn in fovs]
    fov_numbers = [fov_int for fov_int, _ in fovs]  # 1-indexed FOV from group name
    n_cells_per_fov = [int(src[gn]["single_cell_data"].shape[0]) for gn in fov_group_names]
    n_cells = sum(n_cells_per_fov)

    well_group_name = f"{plate_source}__{plate_batch}__{plate_name}__{well}"
    well_group = dst.create_group(well_group_name)
    well_group.attrs["well"] = well
    well_group.attrs["n_fovs"] = n_fovs
    well_group.attrs["n_cells"] = n_cells

    label_shape = (n_fovs, 2, fov_height, fov_width)
    label_arr = well_group.create_array(
        name="label_image",
        shape=label_shape,
        chunks=(1, 1, fov_height, fov_width),
        shards=label_shape,
        dtype=np.uint32,
        compressors=COMPRESSOR,
    )
    for axis0_idx, group_name in enumerate(fov_group_names):
        label_arr[axis0_idx] = np.asarray(src[group_name]["label_image"])

    sample_scd = src[fov_group_names[0]]["single_cell_data"]
    crop_shape = tuple(int(x) for x in sample_scd.shape[1:])
    crop_dtype = sample_scd.dtype
    sample_sci = src[fov_group_names[0]]["single_cell_index"]
    sci_dtype = sample_sci.dtype  # carry source dtype (uint64) verbatim

    if n_cells > 0:
        shard_cells = min(SHARD_CELLS, n_cells)
        scd_arr = well_group.create_array(
            name="single_cell_data",
            shape=(n_cells, *crop_shape),
            chunks=(1, *crop_shape),
            shards=(shard_cells, *crop_shape),
            dtype=crop_dtype,
            compressors=COMPRESSOR,
        )
        # single_cell_index: [src_row, cellid] carried verbatim from broad/
        sci_arr = well_group.create_array(
            name="single_cell_index",
            shape=(n_cells, 2),
            chunks=(n_cells, 2),
            dtype=sci_dtype,
        )
        # fov_index: 1-indexed FOV number per cell (e.g. [1,1,..,2,2,..,9,9,..])
        fov_index_arr = well_group.create_array(
            name="fov_index",
            shape=(n_cells,),
            chunks=(n_cells,),
            dtype=np.uint32,
        )

        offset = 0
        for axis0_idx, group_name in enumerate(fov_group_names):
            n = n_cells_per_fov[axis0_idx]
            if n == 0:
                continue
            scd_arr[offset:offset + n] = np.asarray(src[group_name]["single_cell_data"])
            sci_arr[offset:offset + n] = np.asarray(src[group_name]["single_cell_index"])
            fov_index_arr[offset:offset + n] = fov_numbers[axis0_idx]
            offset += n
        assert offset == n_cells, f"cell offset mismatch: {offset} != {n_cells}"
    else:
        well_group.create_array(
            name="single_cell_data",
            shape=(0, *crop_shape),
            chunks=(1, *crop_shape),
            dtype=crop_dtype,
        )
        well_group.create_array(
            name="single_cell_index",
            shape=(0, 2),
            chunks=(1, 2),
            dtype=sci_dtype,
        )
        well_group.create_array(
            name="fov_index",
            shape=(0,),
            chunks=(1,),
            dtype=np.uint32,
        )


def _validate_round_trip(dst: zarr.Group, n_samples: int) -> None:
    if n_samples <= 0:
        return
    rng = np.random.default_rng(42)
    wells_with_cells = [
        w for w in dst.group_keys() if int(dst[w].attrs["n_cells"]) > 0
    ]
    if not wells_with_cells:
        module_logger.warning("no cells in any well; skipping round-trip validation")
        return
    n_checked = 0
    attempts = 0
    while n_checked < n_samples and attempts < n_samples * 4:
        attempts += 1
        well_group_name = str(rng.choice(wells_with_cells))
        wg = dst[well_group_name]
        n_cells = int(wg.attrs["n_cells"])
        n_fovs = int(wg.attrs["n_fovs"])
        k = int(rng.integers(0, n_cells))
        # cellid is col 1 of single_cell_index (col 0 = src_row, carried verbatim)
        cellid = int(wg["single_cell_index"][k, 1])
        fov_number = int(wg["fov_index"][k])  # 1-indexed
        # Find the axis-0 slot in label_image whose FOV matches.
        # Each FOV occupies a contiguous run in fov_index, so the unique sorted
        # values in fov_index give the axis-0 ordering.
        fov_index_arr = np.asarray(wg["fov_index"])
        unique_fovs_sorted = sorted(set(int(x) for x in fov_index_arr))
        if len(unique_fovs_sorted) != n_fovs:
            raise RuntimeError(
                f"{well_group_name}: n_fovs attr ({n_fovs}) disagrees with "
                f"fov_index unique count ({len(unique_fovs_sorted)})"
            )
        axis0_idx = unique_fovs_sorted.index(fov_number)
        mask = np.asarray(wg["label_image"][axis0_idx, 1]) == cellid
        if not mask.any():
            raise RuntimeError(
                f"round-trip failed: well_group={well_group_name} cell_row={k} "
                f"fov_number={fov_number} (axis0={axis0_idx}) cellid={cellid} "
                f"-> empty mask"
            )
        n_checked += 1
    module_logger.info(f"round-trip OK on {n_checked} sampled cells")


def main(snakemake) -> None:  # noqa: ANN001 — snakemake magic
    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
        datefmt="%d.%m.%y %H:%M:%S",
    )

    src_plate_zarr = Path(snakemake.input.plate_zarr)
    # `plate_v2_dir` is the parent of the zarr (so snakemake's .snakemake_timestamp
    # marker lives next to the zarr rather than inside it). The zarr name matches
    # the plate wildcard, the same convention as Option A.
    dst_plate_zarr = Path(snakemake.output.plate_v2_dir) / f"{snakemake.wildcards.plate}.zarr"

    try:
        reaggregate_plate(src_plate_zarr, dst_plate_zarr)
    except Exception:
        module_logger.exception("reaggregation failed")
        raise

    checkpoint = Path(snakemake.output.checkpoint)
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    checkpoint.write_text(f"reaggregated {dst_plate_zarr}\n")


if __name__ == "__main__":
    if "snakemake" in dir():
        main(snakemake)  # noqa: F821 — injected by snakemake
    else:
        parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
        parser.add_argument("--src-plate-zarr", required=True, type=Path)
        parser.add_argument("--dst-plate-zarr", required=True, type=Path)
        parser.add_argument(
            "--wells",
            default=None,
            help="comma-separated subset of wells (testing only)",
        )
        parser.add_argument("--validate-sample", default=5, type=int)
        parser.add_argument("--log", default=None, help="optional log file path")
        args = parser.parse_args()

        logging.basicConfig(
            filename=args.log,
            level=logging.INFO,
            format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
            datefmt="%d.%m.%y %H:%M:%S",
        )
        if args.log is None:
            logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

        wells_filter = (
            [w.strip() for w in args.wells.split(",")] if args.wells else None
        )
        reaggregate_plate(
            src_plate_zarr=args.src_plate_zarr,
            dst_plate_zarr=args.dst_plate_zarr,
            wells_filter=wells_filter,
            validate_sample=args.validate_sample,
        )
