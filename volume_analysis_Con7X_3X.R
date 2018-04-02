
countData = read.table(file=".../rsem.genes.counts.matrix",header=TRUE,sep="\t",row.names=1,check.names = FALSE)
colData = read.table(file=".../samples_described.txt",head=TRUE,sep="\t",row.names=1,check.names = FALSE)
rownames(colData)=sub("-","_",rownames(colData), fixed=T)
countData <- countData[ , rownames(colData)]
countData = round(countData)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

dds_sub <- dds[, dds$condition %in% c("Con7X", "Con3X")]
dds_sub$condition <- droplevels(dds_sub$condition)
rld_sub <- rlog(dds_sub, blind=FALSE)

plotPCA(rld_sub)
(data <- plotPCA(rld_sub, returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))
library("ggplot2")
shapes <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)
#shapes <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,25)
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())

pdf(".../PCA_7X_3X.pdf")
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
dev.off()


dds_sub$condition <- relevel(dds_sub$condition, ref="Con7X")
dds_sub$condition <- droplevels(dds_sub$condition)
dds_sub <- DESeq(dds_sub)

res05_3Xvs7X <- results(dds_sub, alpha=0.05) 
res05Ordered_3Xvs7X <- res05_3Xvs7X[order(res05_3Xvs7X$padj),]

write.csv(as.data.frame(res05Ordered_3Xvs7X), file=".../Con3Xvs7X_diffexp_adjp005.csv")


### Intersecting diff exp genes with p<0.05 and log2foldchange >1 or <-1 with gene sets ####

data_3Xvs7X <- na.omit(res05Ordered_3Xvs7X)
fltr_3Xvs7X <- data_3Xvs7X[(as.matrix(data_3Xvs7X[,6]) < 0.05)&(as.matrix(data_3Xvs7X[,2]) > 1 | as.matrix(data_3Xvs7X[,2]) < -1),]
dim(fltr_3Xvs7X)

forebrain=read.table(file=".../GO_FOREBRAIN_DEVELOPMENT.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
forebrain_3Xvs7X <- intersect(rownames(fltr_3Xvs7X),rownames(forebrain))
length(forebrain_3Xvs7X)
mt <- assay(rld_sub)[forebrain_3Xvs7X,]
mt <- mt - rowMeans(mt)

library("pheatmap")
library("RColorBrewer")
df <- as.data.frame(colData(dds_sub)[c("condition")])
pheatmap(mt, annotation_col=df, cluster_rows =TRUE, cluster_cols=FALSE, fontsize=6.5,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
pdf(".../forebrain_dev_3Xvs7X_intersect.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =TRUE, cluster_cols=FALSE, fontsize=6.5,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()

kegg_axon_guid=read.table(file=".../kegg_axon_guidance_geneset.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
axonguid_3Xvs7X <- intersect(rownames(fltr_3Xvs7X),rownames(kegg_axon_guid))
length(axonguid_3Xvs7X)
mt <- assay(rld_sub)[axonguid_3Xvs7X,]
mt <- mt - rowMeans(mt)

pheatmap(mt, annotation_col=df, cluster_rows =TRUE, cluster_cols=FALSE, fontsize=6.5,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
pdf(".../kegg_axon_guid_3Xvs7X_intersect.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =TRUE, cluster_cols=FALSE, fontsize=6.5,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()

cent_nervsys_dev=read.table(file=".../GO_CENTRAL_NERVOUS_SYSTEM_DEVELOPMENT.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
cent_nrv_sys_3Xvs7X <- intersect(rownames(fltr_3Xvs7X),rownames(cent_nervsys_dev))
length(cent_nrv_sys_3Xvs7X)
mt <- assay(rld_sub)[cent_nrv_sys_3Xvs7X,]
mt <- mt - rowMeans(mt)

pheatmap(mt, annotation_col=df, cluster_rows =TRUE, cluster_cols=FALSE, fontsize=6.5,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
pdf(".../cent_nervsys_dev_3Xvs7X_intersect.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =TRUE, cluster_cols=FALSE, fontsize=6.5,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()

neuron_dev=read.table(file=".../GO_NEURON_DEVELOPMENT.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
neurondev_3Xvs7X <- intersect(rownames(fltr_3Xvs7X),rownames(neuron_dev))
length(neurondev_3Xvs7X)
mt <- assay(rld_sub)[neurondev_3Xvs7X,]
mt <- mt - rowMeans(mt)

pheatmap(mt, annotation_col=df, cluster_rows =TRUE, cluster_cols=FALSE, fontsize=6.5,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
pdf(".../neuron_dev_3Xvs7X_intersect.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =TRUE, cluster_cols=FALSE, fontsize=6.5,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()

