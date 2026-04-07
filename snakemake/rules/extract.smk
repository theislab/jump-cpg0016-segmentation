rule extract:
    input:
        batch_stack="results/images/{batch}.npy"
    log: "log/extract/{source}/{batch}.log"
    output:
        dir=temp(directory("results/sparcspy/{source}/{batch}/")),
        extraction=temp("results/extraction/{source}/{batch}/extracted_single_cells.h5"),
        segmentation=temp("results/segmentation/{source}/{batch}/input_segmentation.h5")
        # dir=directory("results/sparcspy/{source}/{batch}/"),
        # extraction="results/extraction/{source}/{batch}/extracted_single_cells.h5"
    params:
        config="config/sparcspy.yml"
    priority: 100
    # retries: 5
    resources:
        vram_mb=2000,
        mem_mb=20000,
        tmpdir="results/tmp"
    shell:
        """
        CONDA_PYTHON=/ictstr01/groups/ml01/projects/2023_hackathon23_subcellular_spatial_niklas.schmacke/jump/envs/79ff5f567d72279473b1703565acfd26_/bin/python
        $CONDA_PYTHON scripts/extract.py \
            --batch-stack {input.batch_stack} \
            --project-dir {output.dir} \
            --config {params.config} \
            --metadata {config[samples_meta]} \
            --batch {wildcards.batch} \
            --extraction {output.extraction} \
            --segmentation {output.segmentation} \
            >> {log} 2>&1
        """
