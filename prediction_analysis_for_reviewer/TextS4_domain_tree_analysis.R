
library(ggplot2)
library(reshape2)
library(edgeR)
library(plyr)
library(scales)

library(Rsubread)
library(structSSI)
library(data.table)
library(igraph)

library(cluster)
library(fastcluster)
library(RFLPtools)

library(pheatmap)
library(EDASeq)
library(dynamicTreeCut)
library(moduleColor)
library(knitr)

library(DESeq2)

wd = "./"
colors <- c("#0571b0","#ca0020")

## ------------------------------------------------------------------------
setwd(wd)
# files <- Sys.glob("/home/users/allstaff/tonkin-hill.g/gene_level_var_analysis/alignments_seperate_assemblies_best1/*.sam")
# dom_counts <- featureCounts(files, annot.ext = "./data/NTS_DBL_CIDR_seperate_assembly_domains.SAF"
#                             , isPairedEnd=TRUE, isGTFAnnotationFile=FALSE
#                             , countChimericFragments=TRUE
#                             , allowMultiOverlap=TRUE
#                             , useMetaFeatures=FALSE
#                             , nthreads=20
#                             , countMultiMappingReads=TRUE)
# 
# #combine counts for different runs of the same sample
# colnames(dom_counts$counts) <- gsub(".*_best1\\.","",colnames(dom_counts$counts))
# colnames(dom_counts$counts) <- gsub("_S.*","",colnames(dom_counts$counts))
# temp <- dom_counts
# temp$counts <- t(rowsum(t(temp$counts), group = rownames(t(temp$counts))))
# dim(temp$counts)
# dom_counts <- temp
# # #
# saveRDS(dom_counts, "./data/dom_counts.rds")

dom_counts <- readRDS("./data/dom_counts.rds")

## ------------------------------------------------------------------------
setwd(wd)
dom_counts$counts <- dom_counts$counts[,!(colnames(dom_counts$counts) %in% c('SFD1.CM', 'SFM9', 'SFC.025'))]

## ------------------------------------------------------------------------
setwd(wd)
prev_dom_counts <- data.frame(dom_counts$counts)
prev_dom_counts$domainType <- unlist(lapply(rownames(prev_dom_counts), function(x){
  temp <- unlist(strsplit(x, "_"))
  return(temp[length(temp)-1])
  }))
# prev_dom_counts$domainType <- gsub("\\..*", "", prev_dom_counts$domainType)
prev_dom_counts <- ddply(prev_dom_counts, "domainType", numcolwise(sum))
rownames(prev_dom_counts) <- prev_dom_counts$domainType
prev_dom_counts$domainType <- NULL
prev_dom_counts <- as.matrix(prev_dom_counts)

x_domains <- DGEList(counts=prev_dom_counts)
x_domains <- calcNormFactors(x_domains, method="RLE")
categories <- as.factor(substring(rownames(x_domains$samples),1,1))

## ------------------------------------------------------------------------
plotPCA(cpm(x_domains,normalized.lib.sizes=TRUE, log=T, prior.count=0.5), col=colors[categories], cex=0.8, isLog=TRUE)
title("VAR PCA - Previous Domains")

## ------------------------------------------------------------------------
colData <- data.frame(disease=categories)
rownames(colData) <- colnames(x_domains$counts)
dds_dom <- DESeqDataSetFromMatrix(countData = x_domains$counts
                              , colData = colData
                              , design = ~disease)
dds_dom <- estimateSizeFactors(dds_dom)
dds_dom <- DESeq(dds_dom) 
res <- results(dds_dom,contrast=c("disease","S","I"), pAdjustMethod = "BH")
summary(res, alpha=0.05)

## ------------------------------------------------------------------------
res_domains <- res
resOrdered <- res_domains[order(res_domains$padj),]
resSig <- subset(resOrdered, padj < 0.05)

resSigLFC <- data.frame(subset(resSig, log2FoldChange > 0))

select <- order(rowMeans(counts(dds_dom, normalized=TRUE)), decreasing=TRUE)

nt <- normTransform(dds_dom)
log2.norm.counts <- assay(nt)[select,]
log2.norm.counts <- assay(nt)[rownames(nt) %in% rownames(resSigLFC),]
annotData <- colData

disease        <- colors
names(disease) <- c("non-severe", "severe")
anno_colors <- list(phenotype = disease)
annotData$phenotype[annotData$disease=="I"] <- "non-severe" 
annotData$phenotype[annotData$disease=="S"] <- "severe"
# annotData <- annotData[annotData$disease=="S",]
annotData$disease <- NULL
pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=annotData, annotation_colors = anno_colors)

## ------------------------------------------------------------------------
domain_length_df <- fread("./data/domain_lengths.csv", data.table = FALSE)
domain_length_matrix <- matrix(domain_length_df[,2], ncol=1,
                               dimnames = list(row=domain_length_df[,1], col="len"))
rpkm_dom <- counts(dds_dom, normalized=FALSE)
dom_length_vector <- domain_length_matrix[rownames(rpkm_dom),]
domain_lib_size <- colSums(rpkm_dom)
rpkm_dom <- (t(t(rpkm_dom)/domain_lib_size)/dom_length_vector)*10^9
rpkm_dom <-  rpkm_dom[rownames(rpkm_dom) %in% rownames(resSigLFC),]

# norm_dom_counts <- counts(dds_dom, normalized=TRUE)
# sigUp_norm_counts <- norm_dom_counts[rownames(norm_dom_counts) %in% rownames(resSigLFC),]
# 
# percent_counts <- t(t(sigUp_norm_counts)/colSums(norm_dom_counts))

rpkm_dom <- rpkm_dom+1 #to deal with log transformation
sigUp_norm_counts_df <- melt(rpkm_dom)
colnames(sigUp_norm_counts_df) <- c("Domain Cluster", "Isolate", "RPKM")
sigUp_norm_counts_df$disease[grepl("^S.*",sigUp_norm_counts_df$Isolate)]<-"severe"
sigUp_norm_counts_df$disease[!grepl("^S.*",sigUp_norm_counts_df$Isolate)]<-"non-severe"
colnames(sigUp_norm_counts_df) <- c("Domain Cluster", "Isolate", "RPKM", "Phenotype")

gg <- ggplot(sigUp_norm_counts_df, aes(x=Phenotype, y=`RPKM`
                                       , colour=Phenotype))
gg <- gg + geom_boxplot(outlier.shape = NA, position = "dodge") + geom_jitter(width = 0.2, height = 0)
gg <- gg + facet_wrap(~`Domain Cluster`, nrow=1)
gg <- gg + scale_colour_manual(values = colors) + scale_y_log10()
gg <- gg + theme_bw()
gg

## ------------------------------------------------------------------------
setwd(wd)

select <- order(rowMeans(counts(dds_dom, normalized=TRUE)), decreasing=TRUE)

nt <- normTransform(dds_dom)
log2.norm.counts <- assay(nt)[select,]
log2.norm.counts <- assay(nt)[rownames(nt) %in% c("NTSA", "NTSB"),]

normCounts <- data.frame(t(log2.norm.counts))
normCounts$Sample <- rownames(normCounts)
normCounts$phenotype[as.factor(substring(rownames(normCounts),1,1))=='S'] <-  "severe"
normCounts$phenotype[as.factor(substring(rownames(normCounts),1,1))=='I'] <-  "non-severe"
normCounts <- melt(normCounts, id.vars=c("Sample","phenotype"), variable.name="Type", value.name="Normalised.Count")

gg <- ggplot(normCounts, aes(x=Type, y=Normalised.Count, colour=phenotype))
# gg <- gg + geom_point(position=position_jitter(w=0.1,h=0))
gg <- gg + geom_boxplot()
gg <- gg + theme_bw()
gg <- gg + scale_colour_manual(values = colors)
gg

## ------------------------------------------------------------------------
setwd(wd)
clusters <- read.csv("./data/NTS_DBL_CIDR_seperate_assembly_domains_PreProcessed_hierarchical_clusters.csv")

stopifnot(nrow(dom_counts$counts)==nrow(clusters))

dom_counts_clusters <- merge(dom_counts$counts, clusters, by.x=0 , by.y="Transcript")

## ------------------------------------------------------------------------
setwd(wd)

#first create an edgelist
edges <- list()
edges <- rbind(edges,cbind("root", paste(as.character(clusters[,2]), colnames(clusters)[2], sep="_")))
for (j in 3:(ncol(clusters))){
      edges <- rbind(edges,cbind(paste(as.character(clusters[,j-1]), colnames(clusters)[j-1], sep="_"), paste(as.character(clusters[,j]), colnames(clusters)[j], sep="_")))
  }
edges <- edges[!duplicated(edges),]
edges <- apply(as.matrix(edges), 2,as.character)

#now sum counts at the different clustering heights
results <- data.frame()
pca.df <- data.frame()
all_clusters.matrix <- data.frame()
for (c in colnames(clusters)[2:(ncol(clusters))]) {
  #reduce matrix to cluster at level c
  dom_counts_clusters.dt <- data.table(dom_counts_clusters)
  dom_counts_clusters.dt <- dom_counts_clusters.dt[,lapply(.SD, sum, na.rm=TRUE)
                                                       , by=c
                                                       , .SDcols=colnames(dom_counts_clusters)[2:42]]
  
  dom_counts_clusters.matrix <- as.matrix(dom_counts_clusters.dt)
  dom_counts_clusters.matrix <- dom_counts_clusters.matrix[,2:ncol(dom_counts_clusters.matrix)]
  rownames(dom_counts_clusters.matrix) <- paste(dom_counts_clusters.dt[[1]], c, sep="_")
  all_clusters.matrix <- rbind(all_clusters.matrix, dom_counts_clusters.matrix)
  
  #generate points for PCA
  dgePCA <- DGEList(counts=dom_counts_clusters.matrix)
  colData <- data.frame(disease=as.factor(substring(colnames(dom_counts_clusters.matrix),1,1)))
  rownames(colData) <- rownames(dgePCA$samples)
  ddsPCA <- DESeqDataSetFromMatrix(countData = dgePCA$counts#dom_counts_clusters.matrix
                            , colData = colData
                            , design = ~disease)
  cts <- counts(ddsPCA)
  geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  ddsPCA <- estimateSizeFactors(ddsPCA, geoMeans=geoMeans)
  pca <- plotPCA(varianceStabilizingTransformation(ddsPCA), intgroup=c("disease"), returnData = TRUE)
  pca$level <- c
  pca.df <- rbind(pca.df, pca)
}

## ------------------------------------------------------------------------
setwd(wd)
pca.df$phenotype[pca.df$disease=='S'] <- "severe"
pca.df$phenotype[pca.df$disease=='I'] <- "non-severe"
gg <- ggplot(pca.df, aes(PC1, PC2, color=phenotype)) + scale_color_manual(values=colors)
# gg <- gg + geom_text(aes(label=name))
gg <- gg + geom_point()
gg <- gg + facet_wrap(~level)
gg <- gg + theme_bw() 
gg

## ------------------------------------------------------------------------
setwd(wd)

dom_counts_clusters.dt <- data.table(dom_counts_clusters)
dom_counts_clusters.dt <- dom_counts_clusters.dt[,lapply(.SD, sum, na.rm=TRUE)
                                                       , by='X0.5'
                                                       , .SDcols=colnames(dom_counts_clusters)[2:42]]
  
dom_counts_clusters.matrix <- as.matrix(dom_counts_clusters.dt)
dom_counts_clusters.matrix <- dom_counts_clusters.matrix[,2:ncol(dom_counts_clusters.matrix)]
rownames(dom_counts_clusters.matrix) <- paste(dom_counts_clusters.dt[[1]], c, sep="_")
  
#generate points for PCA
my_domains <- DGEList(counts=dom_counts_clusters.matrix)
my_domains <- calcNormFactors(my_domains, method="RLE")
plotPCA(cpm(my_domains,normalized.lib.sizes=TRUE, log=T, prior.count=0.5), col=colors[categories], cex=0.8, isLog=TRUE)
title("VAR PCA - Domains Clustered at 50% idenitiy")

## ------------------------------------------------------------------------
setwd(wd)

dge_clust <- DGEList(counts=as.matrix(all_clusters.matrix))

#create design matrix
colData <- data.frame(disease=as.factor(substring(rownames(dge_clust$samples),1,1)))
rownames(colData) <- rownames(dge_clust$samples)

cpm(10, mean(dge_clust$samples$lib.size))

keep <- rowSums(cpm(dge_clust)>15) >= 5
dge_clust <- dge_clust[keep,,keep.lib.sizes=FALSE]

## ---- cache=FALSE--------------------------------------------------------
setwd(wd)

#create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = dge_clust$counts
                            , colData = colData
                            , design = ~disease)

#estimate size factors
dds <- estimateSizeFactors(dds)

#run DESeq2 pipeline
dds <- DESeq(dds)
results <- results(dds,contrast=c("disease","S","I"), pAdjustMethod = "BY")
summary(results, alpha=0.05)

## ------------------------------------------------------------------------
setwd(wd)

Sig <- results[!is.na(results$padj),]
Sig <- Sig[Sig$padj<0.05,]

tree <- graph.edgelist(edges)

naNodes <- V(tree)$name[!(V(tree)$name %in% rownames(results))]
na1s = rep(1, length(naNodes))

unadjp <- c(na1s, results$pvalue)
adjp <- c(rep(NA, length(naNodes)), results$padj)
names(unadjp) <- c(naNodes, rownames(results))
names(adjp) <- c(naNodes, rownames(results))


if (!all(names(unadjp) %in% unique(as.vector(edges)))) {
        stop("Names of elements in unadjp do not match names of tree nodes")
    }
if (!all(unique(as.vector(edges)) %in% names(unadjp))) {
        stop("Names of elements in unadjp do not match names of tree nodes")
    }

p.vals <- data.frame(unadjp, adjp = adjp)
rownames(p.vals) <- names(unadjp)

cluster_names <- clusters
for (i in 2:ncol(clusters)){
  cluster_names[,i] <- paste(as.character(cluster_names[,i]), colnames(clusters)[i], sep="_")
}

## ------------------------------------------------------------------------
setwd(wd)

#return the children of a node
get_children <- function(node, edges) {
    if (is.null(node)){
      return(NULL)
    }
    tree.el.tmp <- data.frame(parent = edges[, 1], child = edges[, 2], stringsAsFactors = F)
    children <- tree.el.tmp[which(tree.el.tmp$parent == node), "child"]
    children <- unlist(c(children, lapply(children, get_children, edges)))
    return(children)
}

#return the anscestors of a node
get_parents <- function(node, edges){
  if (is.null(node)){
      return(NULL)
    }
    tree.el.tmp <- data.frame(parent = edges[, 1], child = edges[, 2], stringsAsFactors = F)
    parent <- tree.el.tmp[which(tree.el.tmp$child == node), "parent"]
    parent <- unlist(c(parent, lapply(parent, get_parents, edges)))
    return(parent)
}

#plot the subtree of a node
plot_tree <- function(node, edges, p.vals, return_script=FALSE){
  tree <- graph.edgelist(edges)
  ancestors <- get_parents(node, edges)
  root <- ancestors[length(ancestors)-1] #choose root one level below root of entire tree
  if (ancestors[1]=="root"){
    root <- node
  }
  
  keep <- c(get_children(root, edges), root)
  newTree <- induced.subgraph(tree, vids=V(tree)[V(tree)$name %in% keep])
  p.vals <- p.vals[rownames(p.vals) %in% V(newTree)$name,]
  p.vals <- p.vals[match(unique(as.vector(t(get.edgelist(newTree)))), rownames(p.vals)),]
  
  hyp.tree <- new("hypothesesTree", tree = get.edgelist(newTree),
        p.vals = p.vals, alpha = 0.05)
  
  plot(hyp.tree, adjust = TRUE, return_script=return_script)
}

## ------------------------------------------------------------------------
setwd(wd)

tempSig <- Sig
tempSig$cluster <- as.numeric(gsub(".*_X","",rownames(tempSig)))
tempSig <- tempSig[order(tempSig$cluster),]
SigFilteredResults <- data.frame()
while (nrow(tempSig)>0){
  trow <- tempSig[which.min(tempSig$padj),]
  SigFilteredResults <- rbind(SigFilteredResults, trow)
  t <- rownames(trow)
  tChildren <- get_children(t, edges)
  tParents <- get_parents(t, edges)
  tempSig <- tempSig[!(rownames(tempSig) %in% c(t, tChildren, tParents)),]
  tempSig <- tempSig[order(tempSig$cluster),]
}

resUp <- SigFilteredResults[SigFilteredResults$log2FoldChange>0,]
resDown <- SigFilteredResults[SigFilteredResults$log2FoldChange<0,]

## ------------------------------------------------------------------------
#Collect sequence labels for each cluster
final <- data.frame(SigFilteredResults)
sig_clusters <- as.numeric(gsub("_.*","", rownames(SigFilteredResults)))
levels <- gsub(".*_","", rownames(SigFilteredResults))
significantClusters <- list()
for (i in 1:length(levels)){
  significantClusters[[i]] <- clusters[clusters[,colnames(clusters)==levels[i]]==sig_clusters[i],]$Transcript
}
final$sequences <- I(significantClusters)

getDomain <- function(x){
  x <- x[[1]][1]
  x <- toString(x)
  x <- strsplit(x, "_")
  x <- x[[1]]
  x <- x[length(x)-1]
  return(x)
}
cluster_domain_type <- unlist(lapply(final$sequences, getDomain))
names(cluster_domain_type) <- rownames(final)

## ------------------------------------------------------------------------
# expression <- cpm(dge_clust, normalized.lib.sizes=TRUE, log=TRUE)[rownames(dge_clust) %in% rownames(resUp),]
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)

nt <- normTransform(dds)
expression <- assay(nt)[rownames(nt) %in% rownames(resUp),]

dissim = dist(expression)
dendro = hclust(dissim, method = "complete")

AutoColor = NULL
for (deepSplit in 0:3)
  AutoColor = cbind(AutoColor, labels2colors(cutreeDynamic(dendro = dendro,
                    minClusterSize = 2,
                    cutHeight = NULL, method = "hybrid", deepSplit = deepSplit,
                    distM = as.matrix(dissim))))
AutoLabels = paste("Hybrid 'auto': dS =", c(0:3))

cut <- cutreeDynamic(dendro = dendro,
                    minClusterSize = 2,
                    cutHeight = NULL, method = "hybrid", deepSplit = 0,
                    distM = as.matrix(dissim))

rowLabels <- dendro$labels
# rowLabels <- rowLabels[dendro$order]
rowsAnnot <- data.frame(group=labels2colors(cut))
rownames(rowsAnnot) <- rowLabels


par(mfrow=c(2,1))
par(cex = 1.4);
par(mar = c(0,8.5,2,0));
plot(dendro,labels=F,main="Co-expression clustering with dynamic tree cut", sub="", xlab="")
par(mar = c(1,8.5,0,0));
plotHclustColors(dendro, AutoColor, 
            AutoLabels, 
            main = "")

## ------------------------------------------------------------------------
setwd(wd)

select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)

nt <- normTransform(dds)
log2.norm.counts <- assay(nt)[select,]
log2.norm.counts <- assay(nt)[rownames(nt) %in% rownames(resUp),]
annotData <- colData

disease        <- colors
names(disease) <- c("non-severe", "severe")
anno_colors <- list(phenotype = disease)
annotData$phenotype[annotData$disease=="I"] <- "non-severe"
annotData$phenotype[annotData$disease=="S"] <- "severe"
annotData$disease <- NULL

log2.norm.counts <- log2.norm.counts[dendro$order,]
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=TRUE
         , annotation_col=annotData, annotation_colors = anno_colors
         , clustering_method = "complete", cutree_rows=9
         , labels_row=paste(rownames(log2.norm.counts), cluster_domain_type[rownames(log2.norm.counts)] ,sep=" - ")
         , annotation_row=rowsAnnot, fontsize_row=7)

## ------------------------------------------------------------------------
rpkm_dom_clus <- counts(dds, normalized=FALSE)
rpkm_dom_clus <- rpkm_dom_clus[rownames(rpkm_dom_clus) %in% rownames(resUp),]
dom_length_vector <- domain_length_matrix[cluster_domain_type[rownames(rpkm_dom_clus)],]
rownames(rpkm_dom_clus) <- paste(rownames(rpkm_dom_clus),
                           cluster_domain_type[rownames(rpkm_dom_clus)] ,
                           sep=" - ")
rpkm_dom_clus <- (t(t(rpkm_dom_clus)/domain_lib_size)/dom_length_vector)*10^9

rpkm_dom_clus <- rpkm_dom_clus+1 #to deal with log transformation
sigUp_norm_counts_clus_df <- melt(rpkm_dom_clus)
colnames(sigUp_norm_counts_clus_df) <- c("Domain Cluster", "Isolate", "RPKM")

sigUp_norm_counts_clus_df$disease[grepl("^S.*",sigUp_norm_counts_clus_df$Isolate)]<-"severe"
sigUp_norm_counts_clus_df$disease[!grepl("^S.*",sigUp_norm_counts_clus_df$Isolate)]<-"non-severe"
colnames(sigUp_norm_counts_clus_df) <- c("Domain Cluster", "Isolate", "RPKM", "Phenotype")

combined_rask_new_df <- rbind(sigUp_norm_counts_clus_df, sigUp_norm_counts_df)

gg <- ggplot(combined_rask_new_df, aes(x=Phenotype, y=RPKM
                                       , colour=Phenotype))
gg <- gg + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.05, height = 0)
gg <- gg + facet_wrap(~`Domain Cluster`, nrow=1)
gg <- gg + scale_colour_manual(values = colors) + scale_y_log10()
gg <- gg + theme_bw()
gg <- gg +  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
gg <- gg + theme(panel.spacing.x=unit(0, "lines"),panel.spacing.y=unit(1, "lines"))
gg <- gg + theme(strip.text.x = element_text(angle = 90))
gg

## ------------------------------------------------------------------------
setwd(wd)

nt <- normTransform(dds)
expression <- assay(nt)[rownames(nt) %in% rownames(resDown),]
dissim = dist(expression)
dendro = hclust(dissim, method = "complete")
cut <- cutreeDynamic(dendro = dendro,
                    minClusterSize = 2,
                    cutHeight = NULL, method = "hybrid", deepSplit = 0,
                    distM = as.matrix(dissim))
rowLabels <- dendro$labels
rowsAnnot <- data.frame(group=labels2colors(cut))
rownames(rowsAnnot) <- rowLabels

select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)

nt <- normTransform(dds)
log2.norm.counts <- assay(nt)[select,]
log2.norm.counts <- assay(nt)[rownames(nt) %in% rownames(resDown),]
annotData <- colData

disease        <- colors
names(disease) <- c("non-severe", "severe")
anno_colors <- list(phenotype = disease)
annotData$phenotype[annotData$disease=="I"] <- "non-severe"
annotData$phenotype[annotData$disease=="S"] <- "severe"
annotData$disease <- NULL

log2.norm.counts <- log2.norm.counts[dendro$order,]
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=TRUE
         , annotation_col=annotData, annotation_colors = anno_colors
         , clustering_method = "complete", cutree_rows=10
         , labels_row=paste(rownames(log2.norm.counts), cluster_domain_type[rownames(log2.norm.counts)] ,sep=" - ")
         , annotation_row=rowsAnnot, fontsize_row=6)

## ------------------------------------------------------------------------
setwd(wd)
plotPCA(cpm(dge_clust, normalized.lib.sizes=TRUE, log=T)[rownames(dge_clust) %in% rownames(SigFilteredResults),], col=colors[colData$disease], cex=0.8, isLog=TRUE)

## ------------------------------------------------------------------------
setwd(wd)

#Now print the table
final$padj <- format(final$padj, scientific=TRUE, digits=3)
final$pvalue <- format(final$pvalue, scientific=TRUE, digits=3)
final$sequences <- unlist(lapply(final$sequences, paste, collapse=","))
kable(final, digits=3)

## ------------------------------------------------------------------------
# writeTree <- function(node, directory, edges, p.vals){
#     sink(paste(paste(directory,node,sep=""), "_tree.html", sep=""))
#     cat(plot_tree(node, edges, p.vals, return_script=TRUE))
#     sink()
#     closeAllConnections()
# #     close(paste(paste(directory,node,sep="_"), "_tree.html"))
# }
# plot_tree("780_X0.6", edges, p.vals, return_script=F)
# writeTree("670_X0.6", "/home/users/allstaff/tonkin-hill.g/Dropbox/Var transcriptome manuscript/Figures/VAR differential analysis (4)/domainTrees/", edges, p.vals)
# lapply(rownames(final), writeTree, "/home/users/allstaff/tonkin-hill.g/var_differential_analysis/domain_level/tree_figures/", edges, p.vals)

## ------------------------------------------------------------------------
# plotExpression <- function(cluster, expression, results, outdir){
#   clustExp <- data.frame(expression[rownames(expression)==cluster,])
#   clustExp$Phenotype[as.factor(substring(rownames(clustExp),1,1))=='S'] <-  "severe"
#   clustExp$Phenotype[as.factor(substring(rownames(clustExp),1,1))=='I'] <-  "non-severe"
#   colnames(clustExp) <- c("expression", "phenotype")
#   clustRes <- results[rownames(results)==cluster,]
#   gg <- ggplot(clustExp, aes(x = factor(phenotype), y = expression, fill=factor(phenotype)))
#   gg <- gg + geom_dotplot(binaxis = "y", stackdir = "center", binpositions="all")
#   gg <- gg + theme_bw() + theme(legend.position="none")
#   gg <- gg + scale_fill_manual(values = colors)
#   gg <- gg + labs(title = paste(c("Cluster:", cluster
#                              , "P-value:", format(clustRes$padj, scientific=TRUE, digits=3)
#                              , "LFC:", format(clustRes$log2FoldChange, scientific=TRUE, digits=3))
#                              , sep="", collapse="  ")
#                     , x = "Scaled expression"
#                     , y = "Phenotype"
#                     )
# #   print(gg)
#   ggsave(filename=paste(c(outdir, cluster, "_expression.png"),sep="",collapse=""), plot=gg)
#   }
# 
# 
# 
# lapply(c("91_X0.5", get_parents("91_X0.5", edges), get_children("91_X0.5", edges)), plotExpression, expression, results, "/home/users/allstaff/tonkin-hill.g/var_differential_analysis/domain_level/tree_figures/expression_dot_plots/")

# write.table(results, file="/home/users/allstaff/tonkin-hill.g/var_differential_analysis/domain_level/tree_figures/results.csv", sep=",", quote=F)
# clust <- data.frame(clusters)
# clust$transcript <- clust$Transcript
# clust$Transcript <- NULL
# write.table(clust, file="/home/users/allstaff/tonkin-hill.g/var_differential_analysis/domain_level/tree_figures/segement_clusters_complete.csv", sep=",", quote=F, row.names=F)
#                   
# get_root <- function(node, edges){
#   parents = get_parents(node, edges)
#   if (length(parents)==1){
#     return(node)mai
#   }
#   else{
#     return(parents[length(parents)-1])
#   }
# }
# 
# resSig <- results[!is.na(results$padj),]
# resSig <- resSig[resSig$padj<=0.05,]
# 
# sig <- data.frame(clusters=rownames(SigFilteredResults), root=gsub("_.*", "", unlist(lapply(rownames(SigFilteredResults), get_root, edges))))
# write.table(sig, file="/home/users/allstaff/tonkin-hill.g/var_differential_analysis/domain_level/tree_figures/sig_clusters.csv", sep=",", row.names=FALSE, quote = FALSE)


## ------------------------------------------------------------------------
sessionInfo()

