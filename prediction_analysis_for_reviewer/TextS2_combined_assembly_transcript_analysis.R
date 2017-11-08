library(limma)
library(edgeR)
library(DESeq2)
library(genefilter)
library(EDASeq)
library(RColorBrewer)
library(VennDiagram)
library(knitr)
library(reshape2)

library(pheatmap)
library(dynamicTreeCut)
library(moduleColor)
library(ggplot2)

wd = "./"
colors <- c("#0571b0","#ca0020")

## ------------------------------------------------------------------------
setwd(wd)
load("./data/geneWithoutVivax_and_VARgene_countDB_subread_uH.RData")
var_gene_names <- read.delim("./data/var_gene_names.txt")
var_gene_names <- var_gene_names$X.Gene.ID.
#Load Non-var_genes
fc2 <- fc_genes_aligned_w_Pfhuman_withoutVivax
t <- fc2$counts

#Rename samples to something more sensible
# colnames(t) <- gsub("alignment\\.sam","",colnames(t))
# colnames(t) <- gsub(".*vivax\\.","",colnames(t))
# colnames(t) <- gsub("_S.*","", colnames(t))
nonVar <- t
nonVar.annot <- fc2$annotation

var <- var_counts$counts
var.annot <- var_counts$annotation
# colnames(var) <- gsub(".*alignments\\.","",colnames(var))
# colnames(var) <- gsub("alignment.sam","",colnames(var))
# colnames(var) <- gsub("_L00.*", "", colnames(var))
# colnames(var) <- gsub("_S.*", "", colnames(var))
var_names <- rownames(var)

stopifnot(colnames(var)==colnames(nonVar))
combined_counts <- rbind(nonVar, var)
combined_annot <- rbind(nonVar.annot, var.annot)

## ------------------------------------------------------------------------
combined_counts <- combined_counts[,!(colnames(combined_counts) %in% c('SFD1.CM', 'SFM9', 'SFC.025'))]

## ------------------------------------------------------------------------
stopifnot(rownames(combined_counts)==combined_annot$GeneID)
dge <- DGEList(counts=combined_counts, genes=combined_annot)
dge_combined <- dge

cpm(10, mean(dge_combined$samples$lib.size))

keep <- rowSums(cpm(dge)>1) >= 5
dge <- dge[keep,,keep.lib.sizes=FALSE]

categories <- as.factor(substring(colnames(combined_counts),1,1))

modRUV <- modRUV[match(colnames(dge), rownames(modRUV)),]
stopifnot(colnames(dge)==rownames(modRUV))

## ------------------------------------------------------------------------
categories<- as.factor(substring(rownames(dge$samples),1,1))
dge <- calcNormFactors(dge, method="TMM")
dge$samples$group <- categories
design <- model.matrix(~group, data=dge$samples)

v <- voom(dge, design=design, plot=FALSE)

norm_counts_ring_ruv <- removeBatchEffect(v$E, covariates=modRUV[,3:ncol(modRUV)], design=modRUV[,1:2])
plotPCA(norm_counts_ring_ruv[!(rownames(norm_counts_ring_ruv) %in% c(var_gene_names,var_names)),], col=colors[categories], cex=0.8, isLog=TRUE)
title("Non-VAR PCA")

## ------------------------------------------------------------------------
plotPCA(norm_counts_ring_ruv[rownames(norm_counts_ring_ruv) %in% var_names,], col=colors[categories], cex=0.8, isLog=TRUE)

## ------------------------------------------------------------------------
plotPCA(cpm(dge,normalized.lib.sizes=TRUE, log=T, prior.count=1)[rownames(dge) %in% var_names,], col=colors[categories], cex=0.8, isLog=TRUE)

## ------------------------------------------------------------------------
stopifnot(colnames(dge)==rownames(modRUV))
v2 <- voom(dge, design=modRUV, plot=F)
fit <- lmFit(v2, modRUV)
fit <- eBayes(fit)
deVoom <- decideTests(fit, adjust.method="BH", p.value=0.05)

top_ring_ruv <- topTable(fit, coef=2, p.value=1, sort.by="p", number=Inf, adjust.method="none", confint=TRUE)
top_ring_ruv_var <- top_ring_ruv[rownames(top_ring_ruv) %in% var_names,]
top_ring_ruv_var$adj.P.Val <- p.adjust(top_ring_ruv_var$adj.P.Val, method = "BH")
top_ring_ruv_var <- top_ring_ruv_var[top_ring_ruv_var$adj.P.Val<=0.05,]

voomUp <- rownames(top_ring_ruv_var[top_ring_ruv_var$logFC > 0,])
voomDown <- rownames(top_ring_ruv_var[top_ring_ruv_var$logFC < 0,])
voomUpDown <- c(voomUp, voomDown)

## ------------------------------------------------------------------------
colData <- data.frame(disease=as.factor(substring(rownames(dge_combined$samples),1,1)))
rownames(colData) <- rownames(dge_combined$samples)
dds <- DESeqDataSetFromMatrix(countData = dge_combined$counts[rownames(dge_combined) %in% var_names,]
                              , colData = colData
                              , design = ~disease)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds) 
res <- results(dds,contrast=c("disease","S","I"))

resVar <- res[rownames(res) %in% var_names,]
simpleDE_resVar <- data.frame(resVar)
resOrdered <- resVar
resOrdered <- resOrdered[!is.na(resOrdered$padj),]
resOrdered <- resOrdered[resOrdered$padj<0.05,]

## ------------------------------------------------------------------------
summary(resVar, alpha=0.05)

## ------------------------------------------------------------------------
de_varOnlyGenesUpDown <- rownames(resOrdered)

resUp <- resOrdered[resOrdered$log2FoldChange>0,]
resDown <- resOrdered[resOrdered$log2FoldChange<0,]
DESeqUp <- rownames(resUp)
DESeqDown <- rownames(resDown)
DESeqUpDown <- c(DESeqUp, DESeqDown)

## ------------------------------------------------------------------------
setwd(wd)


select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)

nt <- normTransform(dds)
log2.norm.counts <- assay(nt)[select,]
log2.norm.counts <- assay(nt)[rownames(nt) %in% rownames(resUp),]
# log2.norm.counts <- log2.norm.counts[,colData$disease=="S"]
annotData <- colData

disease        <- colors
names(disease) <- c("non-severe", "severe")
anno_colors <- list(phenotype = disease)
annotData$phenotype[annotData$disease=="I"] <- "non-severe" 
annotData$phenotype[annotData$disease=="S"] <- "severe"
# annotData <- annotData[annotData$disease=="S",]
annotData$disease <- NULL
pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=T, cluster_cols=TRUE, annotation_col=annotData, annotation_colors = anno_colors)

## ------------------------------------------------------------------------
combined_length_df <- fread("./data/combined_assembly_transcript_lengths.csv"
                            , data.table = FALSE)
combined_length_matrix <- matrix(combined_length_df[,2], ncol=1,
                               dimnames = list(row=combined_length_df[,1],
                                               col="len"))
rpkm_combined <- counts(dds, normalized=FALSE)
combined_length_vector <- combined_length_matrix[rownames(rpkm_combined),]
combined_lib_size <- colSums(rpkm_combined)
rpkm_combined <- (t(t(rpkm_combined)/combined_lib_size)/combined_length_vector)*10^9
rpkm_combined <-  rpkm_combined[rownames(rpkm_combined) %in% rownames(resUp),]
rpkm_combined <- rpkm_combined + 1
# combined_rask_new_df <- counts(dds, normalized=TRUE)
# sigUp_norm_counts <- tran_counts[rownames(tran_counts) %in% rownames(resUp),]
# 
# percent_counts <- t(t(sigUp_norm_counts)/colSums(tran_counts))

sigUp_norm_counts_df <- melt(rpkm_combined)
colnames(sigUp_norm_counts_df) <- c("Transcript", "Isolate", "RPKM")
sigUp_norm_counts_df$disease[grepl("^S.*",sigUp_norm_counts_df$Isolate)]<-"severe"
sigUp_norm_counts_df$disease[!grepl("^S.*",sigUp_norm_counts_df$Isolate)]<-"non-severe"
colnames(sigUp_norm_counts_df) <- c("Transcript", "Isolate", "RPKM", "Phenotype")

gg <- ggplot(sigUp_norm_counts_df, aes(x=Phenotype, y=RPKM
                                       , colour=Phenotype))
gg <- gg + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.05, height = 0)
gg <- gg + facet_wrap(~Transcript, nrow=1)
gg <- gg + scale_colour_manual(values = colors) + scale_y_log10()
gg <- gg + theme_bw()
gg <- gg +  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
gg <- gg + theme(panel.spacing.x=unit(0, "lines"),panel.spacing.y=unit(1, "lines"))
gg <- gg + theme(strip.text.x = element_text(angle = 90))
gg

## ------------------------------------------------------------------------
setwd(wd)

select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)

nt <- normTransform(dds)
log2.norm.counts <- assay(nt)[select,]
log2.norm.counts <- assay(nt)[rownames(nt) %in% rownames(resDown),]
# log2.norm.counts <- log2.norm.counts[,colData$disease=="S"]
annotData <- colData

disease        <- colors
names(disease) <- c("non-severe", "severe")
anno_colors <- list(phenotype = disease)
annotData$phenotype[annotData$disease=="I"] <- "non-severe" 
annotData$phenotype[annotData$disease=="S"] <- "severe"
# annotData <- annotData[annotData$disease=="S",]
annotData$disease <- NULL
pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=annotData, annotation_colors = anno_colors)

## ------------------------------------------------------------------------
#Now print the table
resVar$pvalue <- format(resVar$pvalue, scientific=TRUE, digits=3)
resVar$padj <- format(resVar$padj, scientific=TRUE, digits=3)
kable(resVar, digits=3)

## ------------------------------------------------------------------------
sessionInfo()

