rule extract:
    input:
        batch_stack="results/images/{batch}.npy"
    # log: "log/segment/{sample}.log"
    output:
        dir=temp(directory("results/sparcspy/{source}/{batch}/")),
        extraction=temp("results/extraction/{source}/{batch}/extracted_single_cells.h5"),
        segmentation=temp("results/segmentation/{source}/{batch}/input_segmentation.h5")
        # dir=directory("results/sparcspy/{source}/{batch}/"),
        # extraction="results/extraction/{source}/{batch}/extracted_single_cells.h5"
    params:
        config="config/sparcspy.yml",
        debug=False
    priority: 100
    # retries: 5
    resources:
        vram_mb=2000,
        mem_mb=20000,
        tmpdir="results/tmp"
    conda: "../envs/extract.yml"
    script:
        "../scripts/extract.py"
