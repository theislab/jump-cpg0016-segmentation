# Instructions on how to structure the transfer

## Preparation of the transfer per source
During the segmentation, we structured the work by sources. So we now have `n_sources` folders with results on our file system. We will iterate over them. First, we generate a list of "jobs" from a given results folder. Here, a "job" is equal to one fully segmented plate, now stored as a .zarr file. We created a script to automate this generation:

```bash
python 01_extract_joblist_from_results_folder.py \
    --input_path /lustre/groups/ml01/projects/2023_ttreis_segment_JUMP/snakemake/final_source08/snakemake/results \
    --output_path ./source08 \
    --verbosity 1
```

```bash
💡 Valid input path: /ictstr01/groups/ml01/projects/2023_ttreis_segment_JUMP/snakemake/final_source08/snakemake/results
❗ The output path '/ictstr01/home/icb/tim.treis/projects/jump-cpg0016-segmentation/scripts/source08' does not exist.
💡 Creating output path: /ictstr01/home/icb/tim.treis/projects/jump-cpg0016-segmentation/scripts/source08
💡 Valid output path: /ictstr01/home/icb/tim.treis/projects/jump-cpg0016-segmentation/scripts/source08
💡 Dry run: disabled
💡 Verbosity level: 1
💡 Found 1 source(s):
  - source_8
💡 Found 4 batch(es) for source 'source_8'.
💡 Found 56 plate(s) for source 'source_8' in batch 'J3'.
💡 Found 52 plate(s) for source 'source_8' in batch 'J4'.
💡 Found 52 plate(s) for source 'source_8' in batch 'J1'.
💡 Found 56 plate(s) for source 'source_8' in batch 'J2'.
💡 Creating joblist in output directory.
💡 Created 216 job(s).
✅ Processing complete.
```

## Transferring the files
The file transfer will be parallelised across multiple plates using multiple jobs on our HPC, realistically saturating the outgoing traffic. Due to a non-infinite runtime of jobs on our scheduler, these transfer jobs might get cut off when the HPC job dies. So we need to keep track of their status. For this, each job file contains a "status" as it's first line. These can be:
- `todo` -> not yet started
- `ongoing` -> currently being transferred
- `error` -> something went wrong
- `success` -> everything went right
The transfer script also logs its progress in the respective file.

IMPORTANT: We assume that each individual transfer takes less than 12 h, after that the credential for the current transfer would expire. We query a new one after each plate.

```bash
python 04_transfer_plates.py \
    --input_path ./source08 \
    --verbosity 1 
```
```bash
💡 Valid input path: /ictstr01/home/icb/tim.treis/projects/jump-cpg0016-segmentation/scripts/source08
💡 Dry run: disabled
💡 Verbosity level: 1
💡 Found 216 job file(s) in '/ictstr01/home/icb/tim.treis/projects/jump-cpg0016-segmentation/scripts/source08'.
💡 Scanning for 'todo' status...
💡 Found 'todo' job: source_8__J1__A1170383.txt
💡 Uploading files from '/ictstr01/groups/ml01/projects/2023_ttreis_segment_JUMP/snakemake/final_source08/snakemake/results/aggregated/broad/cellpainting-gallery/cpg0016-jump/source_8/workspace/segmentation/cellpose/objects/J1/A1170383' to S3...
<tqdm progress>
```

## Query the progress of the transfer
Because we're tracking the completion of the transfers via static files, we can just analyse them to get a rough overview of our progress:

```bash
python summarise_jobs.py -i ./source08
```
```bash
💡 Valid input path: /ictstr01/home/icb/tim.treis/projects/jump-cpg0016-segmentation/scripts/source08
💡 Dry run: disabled
💡 Verbosity level: 1
💡 Found 216 job file(s).
💡 Total jobs: 216
💡 Todo jobs: 216
💡 Ongoing jobs: 0
✅ Success jobs: 0
❌ Error jobs: 0
```