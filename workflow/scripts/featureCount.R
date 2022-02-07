# This script calculates RPF expression on CDS
library(Rsubread)
library(dplyr)
library(mgsub)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Count RPFs (normalized in RPKM) on CDS for each gene, using `featureCounts`
## run all bams together
samples <- read.table(snakemake@input[["samples"]], header=T)
bamfiles <- paste0("./STAR_align/", as.vector(samples$sample),".bam")

## run one bam file
# bamfiles <- snakemake@input[["bamfile"]]
RPFcounts <- featureCounts(files=bamfiles, annot.ext=snakemake@input[['gtf']],
    isGTFAnnotationFile=TRUE, GTF.featureType="CDS", GTF.attrType="gene_id")

id_length <- RPFcounts$annotation %>% as.data.frame() %>% dplyr::select(GeneID,Length)
rownames(id_length) <- id_length$GeneID
count_table <- merge(id_length, RPFcounts$counts,by="row.names")[,-1]
mapped_reads <- RPFcounts$stat %>% dplyr::filter(Status=="Assigned")
# convert counts to RPKM
values <- mapply('/', count_table %>% summarise(across(starts_with("GSM"), ~./Length*1000*1000000)), mapped_reads[,-1])
rpkm_table <- cbind(count_table[,c(1,2)], values)
colnames(rpkm_table) <- mgsub(colnames(rpkm_table), c("RibosomeProfiling_", ".bam"), c("",""))

write.table(rpkm_table, file=snakemake@output[[1]], quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)