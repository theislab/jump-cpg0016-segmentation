# deactivates as to not take up too much storage
# rule aggregate_sources:
#     input:
#         batches=get_extracted_batches
#     output:
#         h5="results/aggregated/{source}/{source}.h5"
# 
#     log: "log/aggregate/aggregate_{source}.log"
#     conda: "../envs/aggregate_sources.yml"
#     script:
#         "../scripts/aggregate_sources.py"

rule create_broad_structure:
    input:
        extraction_batches=get_extracted_batches,
        segmentation_batches=get_segmentation_batches
    output:
        broad_dir=directory("results/aggregated/broad/cellpainting-gallery/cpg0016-jump/{source}/workspace/segmentation"),
        checkpoint="results/checkpoints/create_broad_structure/{source}.ckpt"
    log: "log/aggregate/broad/create_broad_structure_{source}.log"
    conda: "../envs/aggregate_broad.yml"
    script: "../scripts/create_broad_structure.py"
    

rule aggregate_broad_batch:
    input:
        checkpoint="results/checkpoints/create_broad_structure/{source}.ckpt",
        extraction_path="results/extraction/{source}/{batch}/extracted_single_cells.h5",
        segmentation_path="results/segmentation/{source}/{batch}/input_segmentation.h5",
        broad_dir=directory("results/aggregated/broad/cellpainting-gallery/cpg0016-jump/{source}/workspace/segmentation")
    output:
        checkpoint="results/checkpoints/aggregate_broad_batch/{source}/{batch}.ckpt"
    log: "log/aggregate/broad/aggregate_broad_{batch}_{source}.log"
    conda: "../envs/aggregate_broad.yml"
    script: "../scripts/aggregate_broad_batch.py"

# rule aggregate_broad:
#     input:
#         extraction_batches=get_extracted_batches,
#         segmentation_batches=get_segmentation_batches
#     output:
#         out_dir=directory("results/aggregated/broad/cellpainting-gallery/cpg0016-jump/{source}/workspace/segmentation")

#     log: "log/aggregate/broad/aggregate_broad_{source}.log"
#     conda: "../envs/aggregate_broad.yml"
#     script:
#         "../scripts/aggregate_broad.py"


# Source-level barrier: collapse all per-batch aggregate_broad_batch checkpoints into
# one source-level sentinel. Without this, every reaggregate_to_well job carries N
# per-batch input edges (~159 for source04 × ~277 plates = ~44k stats), which makes
# DAG construction unworkably slow. The barrier rule fires once per source and emits
# a one-line marker file that reaggregate_to_well can depend on as a single edge.
rule aggregate_broad_done:
    input:
        batch_checkpoints=get_aggregate_broad_batch_checkpoints,
    output:
        checkpoint="results/checkpoints/aggregate_broad_done/{source}.ckpt",
    run:
        from pathlib import Path
        Path(output.checkpoint).parent.mkdir(parents=True, exist_ok=True)
        Path(output.checkpoint).write_text(
            f"all aggregate_broad_batch jobs committed for {wildcards.source} "
            f"({len(input.batch_checkpoints)} batches)\n"
        )

# DISABLED 2026-05-20: pipeline pivoted to the per-well aggregated_v2/ deliverable,
# which is produced by `rule reaggregate_to_well` reading directly from broad/.
# This rechunk step (per-FOV, sharded, compressed) is no longer in the chain.
# Code preserved in case we ever need the per-FOV compressed format again.
# rule rechunk_broad:
#     input:
#         batch_checkpoints=get_aggregate_broad_batch_checkpoints,
#         broad_dir=directory("results/aggregated/broad/cellpainting-gallery/cpg0016-jump/{source}/workspace/segmentation")
#     output:
#         compressed_dir=directory("results/aggregated/broad_compressed/cellpainting-gallery/cpg0016-jump/{source}/workspace/segmentation"),
#         checkpoint="results/checkpoints/rechunk_broad/{source}.ckpt"
#     threads: 32
#     log: "log/rechunk/rechunk_{source}.log"
#     conda: "../envs/rechunk.yml"
#     script: "../scripts/rechunk_broad.py"


# Reaggregate per-FOV groups (Option A) to per-well groups (Option B) for one plate.
# Reads from broad/ directly (not broad_compressed/, see disabled rechunk_broad above).
# Writes to sibling tree `aggregated_v2/...` so the source layout is untouched.
# Failure-robust: atomic .tmp -> rename, validation gate, retry on transient failure.
rule reaggregate_to_well:
    input:
        plate_zarr=directory(
            "results/aggregated/broad/cellpainting-gallery/cpg0016-jump/"
            "{source}/workspace/segmentation/{model_dir}/objects/{batch}/{plate}/{plate}.zarr"
        ),
        # Single sentinel for "all per-batch aggregations committed". A plate's
        # FOVs may have been written by multiple snakemake batches into the same
        # broad/{plate}.zarr, so we wait for them all.
        aggregate_done="results/checkpoints/aggregate_broad_done/{source}.ckpt",
    output:
        checkpoint="results/checkpoints/reaggregate_to_well/{source}/{model_dir}/{batch}/{plate}.ckpt",
        # Track the plate-level parent dir, NOT the .zarr itself: snakemake drops
        # a .snakemake_timestamp marker into whichever directory output it's given,
        # and that marker is invalid as a zarr-hierarchy entry. Tracking the parent
        # keeps the marker next to the .zarr and the channel_mapping.json instead.
        plate_v2_dir=directory(
            "results/aggregated_v2/cellpainting-gallery/cpg0016-jump/"
            "{source}/workspace/segmentation/{model_dir}/objects/{batch}/{plate}"
        ),
    log: "log/reaggregate_to_well/{source}/{model_dir}/{batch}/{plate}.log"
    conda: "../envs/aggregate_broad.yml"
    resources:
        mem_mb=4000,
        runtime=120,
    retries: 2
    script: "../scripts/reaggregate_to_well.py"
