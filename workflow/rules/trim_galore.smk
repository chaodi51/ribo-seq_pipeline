rule trim_galore:
    input: 
        get_fastq ## call the function in main Snakefile
    output:
        fq = "trimmed_fq/{sample}_trimmed.fq.gz",
        fq_rep = "trimmed_fq/{sample}.fastq_trimming_report.txt",
    log: "logs/trimmed/{sample}.log"
    threads: 4
    params:
        outdir = trimmed_fq,
        mem = '6G',
        jobName = "trim_galore.{sample}" 
    conda: "../envs/riboseq_env.yaml"
    shell:
        # use BiboToolkit parameters
        "trim_galore --phred33 --length 25 --max_length 34 -e 0.1 -q 30 --stringency 3 --cores {threads} " 
        "--gzip --fastqc -o {params.outdir} {input} &> {log} "