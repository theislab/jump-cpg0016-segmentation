import os
from pathlib import Path
import zarr
import dask.array as da
from tqdm.notebook import tqdm

import joblib
from contextlib import contextmanager

source = "source"
source_num = "08"
cluster_dir = Path("/ictstr01/groups/ml01/projects/2023_ttreis_segment_JUMP/snakemake/")
directory_to_rechunk = cluster_dir / Path(f"final_{source}{source_num}/snakemake/results/aggregated/broad/cellpainting-gallery/cpg0016-jump/source_{int(source_num)}/workspace/segmentation/cellpose/objects")

batch_folders = [f for f in directory_to_rechunk.glob("*") if f.is_dir()]

plate_paths = []
plate_names = []

for batch_folder in batch_folders:
    batch_name = Path(batch_folder.parts[-1])
    batch_folder_contents = [f for f in batch_folder.glob("*") if f.is_dir()]
    plate_names.append([batch_name / Path(f.parts[-1]) for f in batch_folder_contents])

plate_names = [item for sublist in plate_names for item in sublist]

for plate in plate_names:
    plate_id = plate.parts[-1]
    plate_paths.append(directory_to_rechunk / plate / f"{plate_id}.zarr")

# Create a context manager to use tqdm with joblib
@contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar"""
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback

def process_group(store_path, new_store_path, group_name, subgroups):
    """
    Process a single group from the Zarr store.
    """
    src = zarr.open(store_path, mode='r')
    src_group = src[group_name]
    dst_group_path = new_store_path / group_name
    
    for sub in subgroups:
        if sub in src_group:
            src_ds = src_group[sub]
            darray = da.from_zarr(src_ds)
            # Rechunk so that the entire array is one chunk
            darray_single = darray.rechunk(darray.shape)
            dst_ds_path = dst_group_path / sub
            os.makedirs(dst_ds_path, exist_ok=True)
            darray_single.to_zarr(str(dst_ds_path), mode='w', overwrite=True)

subgroups = ['label_image', 'single_cell_data', 'single_cell_index']

n_jobs = os.cpu_count()

for store_path in plate_paths:

    new_store_path = Path(str(store_path).replace('/broad/', '/broad_compressed/'))
    os.makedirs(new_store_path, exist_ok=True)

    src = zarr.open(store_path, mode='r')
    groups = list(src.group_keys())
    
    tasks = [
        joblib.delayed(process_group)(store_path, new_store_path, group_name, subgroups)
        for group_name in groups
    ]
    
    with tqdm(total=len(tasks), desc=f"Processing {store_path.name}") as pbar:
        with tqdm_joblib(pbar):
            joblib.Parallel(n_jobs=n_jobs, verbose=0)(tasks)
