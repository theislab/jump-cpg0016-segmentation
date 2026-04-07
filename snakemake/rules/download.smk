rule download_stacks:
    output: "results/images/{batch}.npy"
    log: "log/download/{batch}.log"
    retries: 5
    resources:
        mem_mb=12000
    shell:
        """
        CONDA_PYTHON=/ictstr01/groups/ml01/projects/2023_hackathon23_subcellular_spatial_niklas.schmacke/jump/envs/507f1bacc8e964624ea0c0fc2f042010_/bin/python
        $CONDA_PYTHON scripts/download_stack.py \
            --batch {wildcards.batch} \
            --metadata {config[samples_meta]} \
            --output {output} \
            --log {log} \
            >> {log} 2>&1
        """
