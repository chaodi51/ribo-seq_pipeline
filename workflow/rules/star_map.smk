# run STAR mapping with only uniquely mapped reads
rule star_map:
    input:
        fq = "cleaned_fq/{sample}.cleaned.fq.gz",
    output:
        "STAR_align/{sample}.bam",
        "STAR_align/{sample}.Log.final.out"
    params:
        genome_dir = STAR_index,
        annotation = gtf_file,
        mem = '10G',
        jobName = "star_map.{sample}" 
    threads: 4
    conda: "../envs/riboseq_env.yaml"
    shell:
    # use the parameters in RiboToolkit
        "STAR --runThreadN {threads} --sjdbGTFfile {params.annotation} --genomeDir {params.genome_dir} " 
        "--outFilterMismatchNmax 2 --quantMode TranscriptomeSAM GeneCounts --outSAMattributes MD NH "
        "--outFilterMultimapNmax 1 --outSAMmapqUnique 255 --outFileNamePrefix STAR_align/{wildcards.sample}. "
        "--outSAMtype BAM SortedByCoordinate --readFilesIn {input} --readFilesCommand gunzip -c && "
        "mv STAR_align/{wildcards.sample}.Aligned.sortedByCoord.out.bam STAR_align/{wildcards.sample}.bam"

rule samtools_index:
    input:
        "STAR_align/{sample}.bam"
    output:
        "STAR_align/{sample}.bam.bai"
    params: 
      mem = '10G',
      jobName = "samtools_index.{sample}"
    conda: "../envs/riboseq_env.yaml"
    shell:
        "samtools index {input} {output}"