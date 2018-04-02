
countData = read.table(file=".../rsem.genes.counts.matrix",header=TRUE,sep="\t",row.names=1,check.names = FALSE)
colData = read.table(file=".../samples_described.txt",head=TRUE,sep="\t",row.names=1,check.names = FALSE)

rownames(colData)=sub("-","_",rownames(colData), fixed=T)

countData <- countData[ , rownames(colData)]
countData = round(countData)

library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds_sub <- dds[, dds$condition %in% c("Con19Z", "Con21Z", "Con22Z")]
dds_sub$condition <- droplevels(dds_sub$condition)
rld_sub <- rlog(dds_sub, blind=FALSE)

### plotting PCA ####

plotPCA(rld_sub)
(data <- plotPCA(rld_sub, returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))
library("ggplot2")
shapes <- c(15,7,17)
#shapes <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,25)
#ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())

ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + scale_shape_manual(values=shapes) + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())

pdf(".../PCA_Con19Z_21Z_22Z.pdf")
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + scale_shape_manual(values=shapes) + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=5) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + scale_shape_manual(values=shapes) + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=6) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + scale_shape_manual(values=shapes) + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
dev.off()


### with Con20Z - 3D only culture of astrocytic cells
dds_sub <- dds[, dds$condition %in% c("Con19Z", "Con21Z", "Con22Z","Con20Z")]
dds_sub$condition <- droplevels(dds_sub$condition)
rld_sub <- rlog(dds_sub, blind=FALSE)

plotPCA(rld_sub)
(data <- plotPCA(rld_sub, returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))
library("ggplot2")
shapes <- c(4,15,7,17)
#shapes <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,25)
#ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())

ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + scale_shape_manual(values=shapes) + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())

pdf(".../with_3Donlycultureofastrocyticcells_PCA_Con19Z_Con20Z_21Z_22Z.pdf")
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + scale_shape_manual(values=shapes) + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=5) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + scale_shape_manual(values=shapes) + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=6) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + scale_shape_manual(values=shapes) + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
dev.off()

selected_genes <- c("GFAP","S100B", "AQP4","AGT","SOX9","VIM","SLC1A3","PAX6","SOX2","HES1","GLI3","SLC17A7","SLC17A6","MAPT","SNAP25")
mat <- assay(rld_sub)[selected_genes,]
library("gplots")
library("RColorBrewer")
heatmap.2(mat,scale="row",trace="none", dendrogram="none", Rowv = FALSE,Colv = FALSE, col = colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255))

pdf(".../astrocytemarkersgenes_profile_with_3Donlycultureofastrocyticcells_PCA_Con19Z_Con20Z_21Z_22Z.pdf")
heatmap.2(mat,scale="row",trace="none", dendrogram="none", Rowv = FALSE,Colv = FALSE, col = colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255))
dev.off()

selected_genes <- c("GFAP","S100B", "AQP4","AGT","SOX9","VIM","SLC1A3","PAX6","SOX2","HES1","GLI3","SLC17A7","SLC17A6","MAPT","SNAP25","CDC20","BIRC5","FAM83D","CDK1")
mat <- assay(rld_sub)[selected_genes,]
heatmap.2(mat,scale="row",trace="none", dendrogram="none", Rowv = FALSE,Colv = FALSE, col = colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255))

pdf(".../astrocyte_&neuronal_&cellcycle_markersgenes_profile_with_3Donlycultureofastrocyticcells_PCA_Con19Z_Con20Z_21Z_22Z.pdf")
heatmap.2(mat,scale="row",trace="none", dendrogram="none", Rowv = FALSE,Colv = FALSE, col = colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255))
dev.off()

### Differential Expression Analysis ###

dds_21Z <- dds[, dds$condition %in% c("Con19Z", "Con21Z")]
dds_21Z$condition <- relevel(dds_21Z$condition, ref="Con19Z")
dds_21Z$condition <- droplevels(dds_21Z$condition)
dds_21Z <- DESeq(dds_21Z)

res05_21Z <- results(dds_21Z, alpha=0.05) 
res05Ordered_21Z <- res05_21Z[order(res05_21Z$padj),]

write.csv(as.data.frame(res05Ordered_21Z), file=".../Con21ZvsCon19Z_diffexp_adjp005.csv")


dds_22Z <- dds[, dds$condition %in% c("Con19Z", "Con22Z")]
dds_22Z$condition <- relevel(dds_22Z$condition, ref="Con19Z")
dds_22Z$condition <- droplevels(dds_22Z$condition)
dds_22Z <- DESeq(dds_22Z)

res05_22Z <- results(dds_22Z, alpha=0.05)
res05Ordered_22Z <- res05_22Z[order(res05_22Z$padj),]

write.csv(as.data.frame(res05Ordered_22Z), file=".../Con22ZvsCon19Z_diffexp_adjp005.csv")

### Checking how many upreg genes compared to Con19Z in Con21Z and Con22Z match ###

data_21Z <- na.omit(res05Ordered_21Z)
upreg_21Zvs19Z_p0.05 <- data_21Z[(as.matrix(data_21Z[,6]) < 0.05)&(as.matrix(data_21Z[,2]) > 0),]
length(rownames(upreg_21Zvs19Z_p0.05))

data_22Z <- na.omit(res05Ordered_22Z)
upreg_22Zvs19Z_p0.05 <- data_22Z[(as.matrix(data_22Z[,6]) < 0.05)&(as.matrix(data_22Z[,2]) > 0),]
length(rownames(upreg_22Zvs19Z_p0.05))

length(intersect(rownames(upreg_21Zvs19Z_p0.05),rownames(upreg_22Zvs19Z_p0.05)))

library(VennDiagram)
grid.newpage()
draw.pairwise.venn(835, 1507, 515, category = c("upreg_21Zvs19Z", "upreg_22Zvs19Z"), fill = c("light blue", "orange"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))

pdf(".../Upreg_genes_venndiagram.pdf")
library(VennDiagram)
grid.newpage()
draw.pairwise.venn(835, 1507, 515, category = c("upreg_21Zvs19Z", "upreg_22Zvs19Z"), fill = c("light blue", "orange"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
dev.off()



