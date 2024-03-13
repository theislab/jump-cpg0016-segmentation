# spatialdata-nik

## How to run the snakemake workflow:

1. Install snakmake via mamba as described [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#full-installation)
2. Install boto3, tqdm, pandas and pandarallel in the environment
    - `mamba activate snakemake`
    - `mamba install boto3 tqdm pandas pandarallel`
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
4. Download metadata for your samples
    - run `python snakemake/scripts/dl.py`
5. Run the snakemake workflow
    - run `snakemake --cores 20 --use-conda --rerun-incomplete --resources gpu=8`
    - adapt number of available cores and gpus to your particular machine
6. The workflow will create a `snakemake/results` directory containing subfolders with the results
    - `results/extraction/<SAMPLE_ID>/single_cells.h5`: File containing the individually extracted cells
    - `results/images/<SAMPLE_ID>.tif`: Stacked tif containing all JUMP channels for a particular sample
    - `results/metadata/`: JUMP Metadata
    - `results/parquet/`: Platewise JUMP Metadata

### More on sample selection

Given the following example `samples.json`:

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

the algorithm would select all samples with the InChIKey `ZGRWVQNYTFGQLL-UHFFFAOYSA-N` and additionally all samples that are both source `source_2` and plate `1053601879`.

### Method to derive `SAMPLE_ID` from JUMP metadata:

The id is created by joining with double underscores. The `Metadata_*` are columns in the metadata tables of JUMP.

```python
f"{row['Metadata_Source']}__{row['Metadata_Batch']}__{row['Metadata_Plate']}__{row['Metadata_Well']}__{row['Metadata_Site']}"
```
