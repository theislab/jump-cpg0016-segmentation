## Segmentation of JUMP cpg0016
The Carpenter-Singh lab at the Broad Institute has recently published the [cpg0016 dataset in which a vast amount of Cell Painting data in response to more than 100k  perturbation was generated](https://github.com/jump-cellpainting/2024_Chandrasekaran_NatureMethods). For current and future downstream ML/DL applications, we have segmented the resulting images and are in the process of uploading them to the Broad's infrastructure. This repository holds the pipeline used for the segmentation and makes it available for inspection and reuse.

## Table of Contents

- [Installation](#installation)
- [Setup](#setup)
- [Running](#running)

## Installation
1. Install snakemake via mamba as described [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#full-installation)
2. Install boto3, tqdm, pandas and pandarallel in the environment
    - `mamba activate snakemake`
    - `mamba install boto3 tqdm pandas pandarallel pyarrow`

## Setup
Running the pipeline is fairly easy once the Python environment is set up since the Cell Painting Gallery is hosted in a way that requires no user accounts, or authentication, and therefore allows the anonymous client of the pipeline to freely download images. When ready, we first have to specify which samples we want to include. This is done with a `samples.json` file that looks as follows:
```json
{
    "samples": [
        {
            "Metadata_InChIKey": "ZGRWVQNYTFGQLL-UHFFFAOYSA-N"
        },
        {
            "Metadata_Source": "source_2",
            "Metadata_Plate": "1053601879"
        }
    ]
}
```
In this example, it would select all samples with the InChIKey `ZGRWVQNYTFGQLL-UHFFFAOYSA-N` and additionally all samples that are both source `source_2` and plate `1053601879`.

The logic to generate these is as follows:
1. Select samples from JUMP that you want to download and process.
    - You can select by the following JUMP metadata columns:
        - Metadata_Source
        - Metadata_Batch
        - Metadata_Plate
        - Metadata_Well
        - Metadata_Site
        - Metadata_InChIKey
        - Metadata_InChI
    - To select samples, modify the `snakemake/config/samples.json` file
        - Each filter in the "samples" list is its own individual filter and resulting samples are added in the end
        - Each metadata column condition in a filter must be fulfilled to return a sample.
        - In `notebooks/generate_example_config.jpynb` we provide examples on how to programmatically generate such a config file.
2. Once specified, we need to download the metadata for the desired samples which is then used during the pipeline.
    - run `python snakemake/scripts/dl.py`

## Running
- Running the pipeline follows standard snakemake logic, for example, using the script in `snakemake/scripts/run_pipeline.sh`. You can optionally specify a directory in which the conda environments for the individual jobs will be created (recommended for debugging purposes).
- Please adapt the number of available cores and GPUs to your particular machine. 

## Further info
The ID is created by joining with double underscores. The `Metadata_*` are columns in the metadata tables of JUMP.

```python
f"{row['Metadata_Source']}__{row['Metadata_Batch']}__{row['Metadata_Plate']}__{row['Metadata_Well']}__{row['Metadata_Site']}"
```


## License
MIT License

[Copyright © 2024, Theislab, Helmholtz Center Munich](./LICENSE)
