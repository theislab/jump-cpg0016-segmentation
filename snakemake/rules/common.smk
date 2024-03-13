import typing
# import subprocess
import pandas as pd

# def run_subprocess(command: str) -> str:
#     try:
#         result = subprocess.run(command.split(), capture_output=True, shell=True, check=True)
#         return result.stdout.decode()
#     except subprocess.CalledProcessError:
#         return "COMMAND_FAILED"
#
# COMMIT_HASH = run_subprocess("git rev-parse HEAD")
# COMMIT_DESCRIPTION = run_subprocess("git describe --all --dirty")
# COMMIT_DIRTY = True if "dirty" in COMMIT_DESCRIPTION else False
# WORKFLOW_HASH = COMMIT_HASH + "-dirty" if COMMIT_DIRTY else COMMIT_HASH

# WORKFLOW_START = datetime.datetime.now(tz=datetime.timezone.utc).strftime('%y%m%d-%H%M%S')

samples = pd.read_parquet(config["samples_meta"])
sources = samples["Metadata_Source"].unique()


def get_agp_url(wildcards): return get_sample_property(wildcards,"s3_OrigAGP")
def get_dna_url(wildcards): return get_sample_property(wildcards,"s3_OrigDNA")
def get_er_url(wildcards): return get_sample_property(wildcards,"s3_OrigER")
def get_mito_url(wildcards): return get_sample_property(wildcards,"s3_OrigMito")
def get_rna_url(wildcards): return get_sample_property(wildcards,"s3_OrigRNA")


def get_batch_images(wildcards) -> typing.List[str]:
    paths = []
    source_samples = samples.loc[
                     (samples["Metadata_Source"] == wildcards.source) &
                     (samples["snakemake_batch"] == wildcards.batch)
    ]
    for _, sample in source_samples.iterrows():
        paths.append(f"results/images/{sample.id}.tif")
    return paths

def get_source_samples(wildcards) -> typing.List[str]:
    source_samples = samples.loc[
        samples["Metadata_Source"] == wildcards.source
    ]
    return source_samples.snakemake_batch.unique().tolist()

def get_extracted_batches(wildcards) -> typing.List[str]:
    paths = [f"results/extraction/{wildcards.source}/{batch}/extracted_single_cells.h5"
             for batch in get_source_samples(wildcards)]
    return paths

def get_segmentation_batches(wildcards) -> typing.List[str]:
    paths = [f"results/segmentation/{wildcards.source}/{batch}/input_segmentation.h5"
             for batch in get_source_samples(wildcards)]
    return paths

def get_sample_property(wildcards, property):
    as_series = samples.loc[samples.id == wildcards.sample][property]
    if len(as_series) == 1:
        return as_series.iloc[0]
    else:
        raise AttributeError(f"Property {property} of {wildcards.sample} did yield multiple results {as_series}")


def get_source_inchi_samples(wildcards) -> typing.List[str]:
    paths = []
    source_samples = samples.loc[samples["Metadata_Source"] == wildcards.source]

    for _, sample in source_samples.iterrows():
        paths.append(f"results/extraction/"
                     f"{sample['Metadata_Source']}/"
                     f"{sample['Metadata_InChIKey']}/"
                     f"{sample['id']}/"
                     f"{sample['id']}.h5")

    return paths

def get_source_inchi_status(wildcards) -> typing.List[str]:
    paths = []
    source_samples = samples.loc[samples["Metadata_Source"] == wildcards.source]

    for _, sample in source_samples.iterrows():
        paths.append(f"results/extraction/"
                     f"{sample['Metadata_Source']}/"
                     f"{sample['Metadata_InChIKey']}/"
                     f"{sample['id']}/"
                     f"segmentation_status.txt")

    return paths

# def get_run_stamp(wildcards) -> str:
#     return WORKFLOW_HASH
