library(riboWaltz)
library(ggplot2)
library(ggpubr)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

annotation_dt <- create_annotation(gtfpath=snakemake@input[["gtf"]], dataSource="UCSC", organism="Mus musculus")

# load data
reads_list <- bamtolist(bamfolder=snakemake@params[["bam_folder"]], annotation=annotation_dt)
# filter by read length
# filtered_list <- length_filter(data = reads_list, length_filter_mode = "custom", length_filter_vector = 26:32)

# P-site offset calculation
POs <- psite(reads_list, flanking=6, extremity="auto")

# update reads_list with p-site info
reads_psite_list <- psite_info(reads_list, POs)

# ------ Plots to overview the data ------- #
## 1.RPF length distribution
all_samples = names(reads_list)
figs=list()
for(sample_name in all_samples){
  length_dist <- rlength_distr(reads_list, sample=sample_name)
  figs <- c(figs, list(length_dist[["plot"]]))
}
p <- ggarrange(plotlist=figs, nrow=2, ncol=2)
pdf(snakemake@output[["RPF_length_plot"]],14,12)
print(p)
dev.off()

## 2.metaheatmaps displays the abundance of the 5' and 3' end of reads on CDSs
figs=list()
for(sample_name in all_samples){
  ends_heatmap <- rends_heat(reads_list, annotation_dt, sample=sample_name, utr5l = 25, cdsl = 50, utr3l = 25)
  figs <- c(figs, list(ends_heatmap[["plot"]] + theme(plot.title=element_text(size=20))))
}
p <- ggarrange(plotlist=figs, nrow=2, ncol=2)
pdf(snakemake@output[["RPF_ends_heatmap_plot"]],20,10)
print(p)
dev.off()

## 3.P-sites per region (5’UTRs, CDSs and 3’UTRs)
figs=list()
for(sample_name in all_samples){
  psite_region <- region_psite(reads_psite_list, annotation_dt, sample=sample_name)
  figs <- c(figs, list(psite_region[["plot"]]))
}
p <- ggarrange(plotlist=figs, nrow=2, ncol=2)
pdf(snakemake@output[["psite_region_plot"]],14,12)
print(p)
dev.off()

## 4.P-sites/RPF signal on three reading frames for 5’UTRs, CDSs and 3’UTRs
### 4.1 by length
figs=list()
for(sample_name in all_samples){
  frames_stratified <- frame_psite_length(reads_psite_list, sample=sample_name, region = "all")
  figs <- c(figs, list(frames_stratified[["plot"]]))
}
p <- ggarrange(plotlist=figs, nrow=2, ncol=2)
pdf(snakemake@output[["Psite_signal_bylength_inframes_plot"]],14,12)
print(p)
dev.off()
### 4.2 in total
figs=list()
for(sample_name in all_samples){
  frames <- frame_psite(reads_psite_list, sample =sample_name, region = "all")
  figs <- c(figs, list(frames[["plot"]]))
}
p <- ggarrange(plotlist=figs, nrow=2, ncol=2)
pdf(snakemake@output[["Psite_signal_total_inframes_plot"]],14,12)
print(p)
dev.off()

## 5.metaplot to show trinucleotide periodicity along CDSs
figs=list()
for(sample_name in all_samples){
  metaprofile <- metaprofile_psite(reads_psite_list, annotation_dt, sample=sample_name, utr5l=20, cdsl=50, utr3l=20,    plot_title="sample.transcript")
  figs <- c(figs, list(metaprofile[[paste0("plot_",sample_name)]] + theme(plot.title=element_text(size=20))))
}
p <- ggarrange(plotlist=figs, nrow=2, ncol=2)
pdf(snakemake@output[["trinucleotide_periodicity_metaprofile_plot"]],24,12)
print(p)
dev.off()

## 6. codon usage
figs=list()
for(sample_name in all_samples){
  cu_barplot <- codon_usage_psite(reads_psite_list, annotation_dt, sample = sample_name, fastapath = snakemake@input[["transcript_fasta"]], fasta_genome = FALSE, frequency_normalization = TRUE) 
  figs <- c(figs, list(cu_barplot[["plot"]] + ggtitle(sample_name) +  theme(plot.title=element_text(size=20))))
}
p <- ggarrange(plotlist=figs, nrow=2, ncol=2)
pdf(snakemake@output[["codon_usage_plot"]],24,12)
print(p)
dev.off()


