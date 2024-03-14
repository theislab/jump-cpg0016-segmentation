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


rule aggregate_broad:
    input:
        extraction_batches=get_extracted_batches,
        segmentation_batches=get_segmentation_batches
    output:
        out_dir=directory("results/aggregated/broad/cellpainting-gallery/cpg0016-jump/{source}/workspace/segmentation")

    log: "log/aggregate/broad/aggregate_broad_{source}.log"
    conda: "../envs/aggregate_broad.yml"
    script:
        "../scripts/aggregate_broad.py"
