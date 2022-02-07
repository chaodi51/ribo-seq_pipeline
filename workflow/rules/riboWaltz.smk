rule samtools_sort:
    input:
        "STAR_align/{sample}.Aligned.toTranscriptome.out.bam"
    output:
        "STAR_align/{sample}.Aligned.toTranscriptome.out.sorted.bam"
    params: 
      mem = '20G',
      jobName = "samtools_sort.{sample}"
    conda: "../envs/riboseq_env.yaml"
    shell:
        "samtools sort {input} -o {output}"

## symbolic link of bam files to shorter names, b/c riboWaltz plots directly use the sample names
rule link_bam:
    input:
        expand("STAR_align/{sample}.Aligned.toTranscriptome.out.sorted.bam", sample=SAMPLES)
    # output:
        # directory("STAR_align_transcriptome/") ## must with / in this case
    params:
        mem = '1G',
        jobName = "link_bam",
        outdir = "STAR_align_transcriptome"
    log: 
        "logs/link_bam.txt"
    threads: 1
    shell:
        '''
        rm -f {log} 
        for i in {input}; do
            sampleName="$(basename $i .Aligned.toTranscriptome.out.sorted.bam)" 2>> {log};
            sample_shortName=`echo $sampleName | sed 's/RibosomeProfiling_//g'` 2>> {log};
            ln -sf ../$i {params.outdir}/$sample_shortName.bam 2>> {log}
        done
        '''

rule riboWaltz:
    input:
        gtf = gtf_file,
        transcript_fasta = transcript_fa
    output:
        RPF_length_plot = "../results/riboWaltz/RPF_length.pdf",
        RPF_ends_heatmap_plot = "../results/riboWaltz/RPF_ends_heatmap.pdf",
        psite_region_plot = "../results/riboWaltz/psite_region.pdf",
        Psite_signal_bylength_inframes_plot = "../results/riboWaltz/Psite_signal_bylength_inframes.pdf",
        Psite_signal_total_inframes_plot = "../results/riboWaltz/Psite_signal_total_inframes.pdf",
        trinucleotide_periodicity_metaprofile_plot = "../results/riboWaltz/trinucleotide_periodicity_metaprofile.pdf",
        codon_usage_plot = "../results/riboWaltz/codon_usage.pdf"
    params:
        bam_folder = "STAR_align_transcriptome/",
        mem = '20G',
        jobName = 'run_riboWaltz'
    # conda:
    #     "../envs/riboWaltz.yaml"
    log:
        "logs/riboWaltz.log"
    threads: 1
    script:
        "../scripts/riboWaltz.R"