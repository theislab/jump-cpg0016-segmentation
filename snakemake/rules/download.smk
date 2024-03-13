rule download_stacks:
    output: "results/images/{batch}.npy"
    log: "log/download/{batch}.log"
    conda: "../envs/download.yml"
    retries: 5
    resources:
        mem_mb=12000
    script:
        "../scripts/download_stack.py"
