## remove the contaminations of rRNAs, tRNA and snRNA to using Bowtie
## with a maximum of two mismatchs (-v 2) by default. tRNA sequences were downloaded from the GtRNAdb database. 
## rRNA and snRNA sequences were retrieved from noncoding RNA annotations in Ensembl

rule bowtie_ncRNAs:
    input:
        fq = "trimmed_fq/{sample}_trimmed.fq.gz",
    output:
        "cleaned_fq/{sample}.cleaned.fq.gz"
    log: "logs/bowtie_ncRNA/{sample}.log"
    params:
        index = bowtie_index_ncRNAs,
        mem = '10G',
        jobName = "bowtie_ncRNAs.{sample}" 
    threads: 4
    conda: "../envs/riboseq_env.yaml"
    shell:
        '''
        bowtie -x {params.index} --threads {threads} {input} --un cleaned_fq/{wildcards.sample}.cleaned.fq -S cleaned_fq/{wildcards.sample}.sam &> {log}
        gzip cleaned_fq/{wildcards.sample}.cleaned.fq 
        rm -f cleaned_fq/{wildcards.sample}.sam
        '''
        
rule bowtie_ncRNA_stat:       
    input: 
        expand("logs/bowtie_ncRNA/{sample}.log", sample=SAMPLES)
    output:  
        "report/mapped_reads_to_ncRNA.txt"
    params:
        mem = '1G',
        jobName = "bowtie_ncRNA_stat"
    shell:
        '''
        rm -f {output};
        for i in {input}; do
            sampleName="$(basename $i .log)";
            cat $i | grep "reads processed" | awk -F": " '{{print $2}}'  > foo1;
            cat $i | grep "reads with at least one alignment" | awk -F": " '{{print $2}}' | sed 's/ /\t/g;s/(//g;s/%)//g' > foo2;
            paste foo1 foo2 | awk '{{print "'$sampleName'\t"$1"\t"$2"\t"$3}}' >> {output};
        done
        sed -i '1isample\ttotal_reads\tmapped_reads\t%mapped' {output};
        rm -f foo*
        '''

