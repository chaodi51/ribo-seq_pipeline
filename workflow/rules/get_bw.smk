## convert the bam file to bigwig format (signal normalized to RPM) for viewing the data on Genome Browse
rule get_bw:
    input: 
        bam = "STAR_align/{sample}.bam",
        rpm_factors = 'report/rpm_factor.txt'
    output: 
        "bw_rpm/{sample}.bw"
    log: "log/get_bw/{sample}.log"
    threads: 1
    params:
        mem = '10G',
        # rpmFactor = lambda wildcards: sample_factor_dict[wildcards.sample],
        jobName = "get_bw.{sample}"
    conda: "../envs/rnaseq_env.yaml"
    shell:
        '''
        factor=`cat {input.rpm_factors} |  grep {wildcards.sample} | cut -f2`;
        genomeCoverageBed -split -bg -ibam {input.bam} -scale $factor 1> bw_rpm/{wildcards.sample}.bg 2> {log};
        bedtools sort -i bw_rpm/{wildcards.sample}.bg 1> bw_rpm/{wildcards.sample}.sort.bg 2>> {log};
        bedGraphToBigWig bw_rpm/{wildcards.sample}.sort.bg STAR_index/chrNameLength.txt {output} 2>> {log};
        rm bw_rpm/{wildcards.sample}.bg bw_rpm/{wildcards.sample}.sort.bg
        '''

