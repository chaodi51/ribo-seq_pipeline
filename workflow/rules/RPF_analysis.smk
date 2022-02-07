rule RPF_exp:
    input:
        # bamfile = "STAR_align/{sample}.bam",
        samples = "sample_contrast.tsv",
        gtf = gtf_file
    output:
        "../results/RPF.cds_exp.tsv"
    log:
        "logs/RPF_exp.log"   
    threads: 1     
    params: 
      mem = '20G',
      jobName = "RPF_exp"
    script:
        "../scripts/featureCount.R"


rule ribo_pause:
    input:
        bamfile = "STAR_align/{sample}.Aligned.toTranscriptome.out.sorted.bam",
        transcript_fasta = transcript_fa 
    output:
        "../results/ribo_pause/{sample}.pause_sites.tsv"
    threads: 1
    params: 
      mem = '20G',
      jobName = "ribo_pause.{sample}",
    #conda: "../envs/riboseq_env.yaml"
    shell:
        '''
        # module load perl/5.26.1,  then need to export the perl path, or use the full perl path directly
        mkdir -p ../results/tmp
        /cm/shared/apps_chop/perl/5.26.1/bin/perl ~/public/tools/Pausepred_offline/offline_pausepred.pl \
        {input.bamfile} \
        1000 20 {input.transcript_fasta} \
        28,29,30 10 50 50 0,0,0 ../results/tmp/{wildcards.sample}.pause_site.txt  
        sed 's/,/\t/g' ../results/tmp/{wildcards.sample}.pause_site.txt > {output}

        '''

# use RiboCode, prepare transcript annotation in ~public/genome/UCSC/mm10
rule activeORF:
    input:
        bamfile = "STAR_align/{sample}.Aligned.toTranscriptome.out.sorted.bam",
        genome_bamfile = "STAR_align/{sample}.bam",
        genome = genome_fa,
        gtf = "/home/dic/public/genomes/UCSC/mm10/mm10.refGene.updated_riboCode.gtf"
    output:
        config = "../results/RiboCode/{sample}/metaplots_pre_config.txt",
        orf_tab = "../results/RiboCode/{sample}.txt",
        orf_count = "../results/RiboCode/{sample}.orf_count"
    log:
        "logs/activeORF/{sample}.log"
    threads: 1     
    params: 
      mem = '20G',
      jobName = "activeORF.{sample}",
      annot_dir = "/home/dic/public/genomes/UCSC/mm10/RiboCode_annot",
      metaplot = "../results/RiboCode/{sample}/metaplots",
      outname = "../results/RiboCode/{sample}",
    conda: "../envs/riboseq_env.yaml"  
    shell:
        '''
        rm -f {log}
        # (1). Preparing the transcripts annotation files: aleady done by 
        # /home/dic/public/genomes/UCSC/mm10/RiboCode_prep.sh
        # prepare_transcripts -g <gencode.v19.annotation.gtf> -f <hg19_genome.fa> -o <RiboCode_annot>
        
        #####
        # (2). Selecting the length range of the RPF reads and identify the P-site locations:
        mkdir -p {params.outname}
        metaplots -a {params.annot_dir} -r {input.bamfile} -s yes -o {params.metaplot} &> {log}
        
        ####
        # (3). Detecting translated ORFs using the ribosome-profiling data:
        RiboCode -a {params.annot_dir} -c {output.config} -l no -g -o {params.outname} &>> {log}
        # (4). (optional) Plotting the P-sites densities of predicted ORFs
        #plot_orf_density -a <RiboCode_annot> -c <config.txt> -t (transcript_id) -s (ORF_gstart) -e (ORF_gstop) 

        ####
        # (5). (optional) Counting the number of RPF reads aligned to ORFs
        ORFcount -g {params.outname}.gtf -r {input.genome_bamfile} -f 15 -l 5 -e 100 -m 26 -M 34 -o {output.orf_count} &>> {log}
        '''
