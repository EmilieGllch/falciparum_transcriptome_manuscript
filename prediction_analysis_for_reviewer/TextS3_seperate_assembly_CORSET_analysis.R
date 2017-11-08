
library(DESeq2)
library(ggplot2)
library(edgeR)
library(knitr)
library(pheatmap)
library(EDASeq)

wd <- "./"
colors <- c("#0571b0","#ca0020")

## ------------------------------------------------------------------------
setwd(wd)

clusters <- read.delim("./data/clusters-0.3.txt", header=FALSE)
colnames(clusters) <- c("Transcript", "Cluster")

transcriptCounts <- read.table("./data/counts-0.3.txt", header=TRUE, sep="\t", row.names=1)
colnames(transcriptCounts) <- gsub("_S.*","",colnames(transcriptCounts))

## ------------------------------------------------------------------------
setwd(wd)
transcriptCounts <- transcriptCounts[,!(colnames(transcriptCounts) %in% c('SFD1.CM', 'SFM9', 'SFC.025'))]

## ------------------------------------------------------------------------
dge <- DGEList(counts=transcriptCounts)
cpm(10, mean(dge$samples$lib.size))
keep <- rowSums(cpm(dge)>180) >= 5
dge <- dge[keep,,keep.lib.sizes=FALSE]
categories <- as.factor(substring(rownames(dge$samples),1,1))
plotPCA(cpm(calcNormFactors(dge, method="RLE"), normalized.lib.sizes=TRUE, log=T), col=colors[categories], cex=0.8, isLog=TRUE)

## ------------------------------------------------------------------------
#create design matrix
colData <- data.frame(disease=as.factor(substring(colnames(transcriptCounts),1,1)))
rownames(colData) <- colnames(transcriptCounts)

#create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = transcriptCounts
                            , colData = colData
                            , design = ~disease)

##estimate Sizefactors
dds <- estimateSizeFactors(dds)

#run DESeq2 pipeline
dds <- DESeq(dds)
res <- results(dds,contrast=c("disease","S","I"), pAdjustMethod = "BH")
summary(res, alpha=0.05)

Sig <- res[!is.na(res$padj),]
Sig <- Sig[Sig$padj<0.05,]

SigUp <- Sig[Sig$log2FoldChange>0,]
SigDown <- Sig[Sig$log2FoldChange<0,]

## ------------------------------------------------------------------------
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)

nt <- normTransform(dds)
log2.norm.counts <- assay(nt)[select,]
log2.norm.counts <- assay(nt)[rownames(nt) %in% rownames(SigUp),]

disease <- colors
names(disease) <- c("non-severe", "severe")
anno_colors <- list(phenotype = disease)
colData$phenotype[colData$disease=="I"] <- "non-severe"
colData$phenotype[colData$disease=="S"] <- "severe"
colData$disease <- NULL

pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=colData, annotation_colors = anno_colors)

## ------------------------------------------------------------------------
clustered_length_df <- fread("./data/sig_clusters_combined_length.csv"
                            , data.table = FALSE)
clustered_length_matrix <- matrix(clustered_length_df[,2], ncol=1,
                               dimnames = list(row=clustered_length_df[,1],
                                               col="len"))
rpkm_clustered <- counts(dds, normalized=FALSE)
clustered_lib_size <- colSums(rpkm_clustered)

rpkm_clustered <- rpkm_clustered[rownames(rpkm_clustered) %in% rownames(SigUp),]

clustered_length_vector <- clustered_length_matrix[rownames(rpkm_clustered),]

rpkm_clustered <- (t(t(rpkm_clustered)/clustered_lib_size)/clustered_length_vector)*10^9
rpkm_clustered <- rpkm_clustered + 1

# tran_counts <- counts(dds, normalized=TRUE)
# sigUp_norm_counts <- tran_counts[rownames(tran_counts) %in% rownames(SigUp),]
# 
# percent_counts <- t(t(sigUp_norm_counts)/colSums(tran_counts))


sigUp_norm_counts_df <- melt(rpkm_clustered)
colnames(sigUp_norm_counts_df) <- c("Transcript", "Isolate", "RPKM")
sigUp_norm_counts_df$disease[grepl("^S.*",sigUp_norm_counts_df$Isolate)]<-"severe"
sigUp_norm_counts_df$disease[!grepl("^S.*",sigUp_norm_counts_df$Isolate)]<-"non-severe"
colnames(sigUp_norm_counts_df) <- c("Transcript Clusters", "Isolate", "RPKM", "Phenotype")

gg <- ggplot(sigUp_norm_counts_df, aes(x=Phenotype, y=RPKM
                                       , colour=Phenotype))
gg <- gg + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.05, height = 0)
gg <- gg + facet_wrap(~`Transcript Clusters`, nrow=1)
gg <- gg + scale_colour_manual(values = colors) + scale_y_log10()
gg <- gg + theme_bw()
gg <- gg +  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
gg <- gg + theme(panel.spacing.x=unit(0, "lines"),panel.spacing.y=unit(1, "lines"))
gg <- gg + theme(strip.text.x = element_text(angle = 90))
gg

## ------------------------------------------------------------------------
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)

nt <- normTransform(dds)
log2.norm.counts <- assay(nt)[select,]
log2.norm.counts <- assay(nt)[rownames(nt) %in% rownames(SigDown),]


pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=colData, annotation_colors = anno_colors)

## ------------------------------------------------------------------------
#Collect sequence labels for each cluster
getTrans <- function(x) {
  return(paste(clusters$Transcript[clusters$Cluster==x], collapse=","))
}

final <- res
final$transcripts <- lapply(rownames(final), getTrans)
#Now print the table
final <- final[order(final$padj),]
final$padj <- format(final$padj, scientific=TRUE, digits=3)
final$pvalue <- format(final$pvalue, scientific=TRUE, digits=3)
kable(final, digits=3)

## ------------------------------------------------------------------------
sessionInfo()

