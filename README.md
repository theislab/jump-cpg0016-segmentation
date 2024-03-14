## Segmentation of JUMP cpg0016
The Carpenter-Singh lab at the Broad institute has recently published the [cpg0016 dataset in which a vast amount of Cell Painting data in response to more than 100k  perturbation was generated](https://github.com/jump-cellpainting/2024_Chandrasekaran_NatureMethods). For current and future downstream ML/DL applications, we have segmented the resulting images and are in the process of uploading them to the Broad's infrastrucutre. This repository holds the pipeline used for the segmentation and makes it available for inspection and reuse.

## Table of Contents

- [Installation](#install)
- [Usage](#usage)
- [License](#license)

## Installation
1. Install snakemake via mamba as described [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#full-installation)
2. Install boto3, tqdm, pandas and pandarallel in the environment
    - `mamba activate snakemake`
    - `mamba install boto3 tqdm pandas pandarallel pyarrow`
3. Select samples from JUMP that you want to download and process.
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
4. Download metadata for your samples
    - run `python snakemake/scripts/dl.py`
5. Run the snakemake workflow
    - run the pipeline, for example using the script in `snakemake/scripts/run_pipeline.sh`. You can optionally specify a directory in which the conda environments for the individual jobs will be created (recommended for debugging purposes).
    - Please adapt number of available cores and GPUs to your particular machine

## License
MIT License

[Copyright © 2024, Theislab, Helmholtz Center Munich](./LICENSE)
