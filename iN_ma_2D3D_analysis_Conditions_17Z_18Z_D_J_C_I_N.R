
countData = read.table(file=".../rsem.genes.counts.matrix",header=TRUE,sep="\t",row.names=1,check.names = FALSE)
colData = read.table(file=".../samples_described.txt",head=TRUE,sep="\t",row.names=1,check.names = FALSE)

rownames(colData)=sub("-","_",rownames(colData), fixed=T)
countData <- countData[ , rownames(colData)]
countData = round(countData)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)


#dds_sub <- dds[, dds$condition %in% c("ConJ", "Con14X", "Con4X", "Con1X", "ConT")]
#dds_sub <- dds[, dds$condition %in% c("Con18Z", "ConK", "Con2X")]
#dds_sub <- dds[, dds$condition %in% c("Con18Z", "ConK")]

dds_sub <- dds[, dds$condition %in% c("Con17Z", "Con18Z", "ConD", "ConJ")]
dds_sub$condition <- droplevels(dds_sub$condition)
rld_sub <- rlog(dds_sub, blind=FALSE)

plotPCA(rld_sub)
(data <- plotPCA(rld_sub, returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))

library("ggplot2")
# shapes <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)
# ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())

shapes <- c(9,13,8,11)
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + scale_shape_manual(values=shapes) + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())

pdf(".../PCA_Con_17Z_18Z_D_J_manualshape_size4.pdf")
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + scale_shape_manual(values=shapes) + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
dev.off()

pdf(".../PCA_Con_17Z_18Z_D_J_manualshape_size4.5.pdf")
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4.5) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + scale_shape_manual(values=shapes) + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
dev.off()

pdf(".../PCA_Con_17Z_18Z_D_J_manualshape_size5.pdf")
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=5) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + scale_shape_manual(values=shapes) + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
dev.off()

pdf(".../PCA_Con_17Z_18Z_D_J_manualshape_size6.pdf")
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=6) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + scale_shape_manual(values=shapes) + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
dev.off()


### Differential Expression Analysis ###

dds_ConJvs18Z <- dds[, dds$condition %in% c("Con18Z", "ConJ")]
dds_ConJvs18Z$condition <- relevel(dds_ConJvs18Z$condition, ref="Con18Z")
dds_ConJvs18Z$condition <- droplevels(dds_ConJvs18Z$condition)
dds_ConJvs18Z <- DESeq(dds_ConJvs18Z)

res05_ConJvs18Z <- results(dds_ConJvs18Z, alpha=0.05) 
res05Ordered_ConJvs18Z <- res05_ConJvs18Z[order(res05_ConJvs18Z$padj),]

write.csv(as.data.frame(res05Ordered_ConJvs18Z), file=".../ConJvsCon18Z_diffexp_adjp005.csv")


### 

dds_ConIvsN <- dds[, dds$condition %in% c("ConN", "ConI")]

dds_ConIvsN$condition <- relevel(dds_ConIvsN$condition, ref="ConN")
dds_ConIvsN$condition <- droplevels(dds_ConIvsN$condition)
dds_ConIvsN <- DESeq(dds_ConIvsN)

res05_ConIvsN <- results(dds_ConIvsN, alpha=0.05) 
res05Ordered_ConIvsN <- res05_ConIvsN[order(res05_ConIvsN$padj),]


## upregulated genes in ConJ for ConJvs18Z (with adjp<0.05)

data_ConJvs18Z <- na.omit(res05Ordered_ConJvs18Z)
upreg_ConJvs18Z_p0.05 <- data_ConJvs18Z[(as.matrix(data_ConJvs18Z[,6]) < 0.05)&(as.matrix(data_ConJvs18Z[,2]) > 0),]
dim(upreg_ConJvs18Z_p0.05)

## upregulated genes in ConI for ConIvsN (with adjp<0.05)

data_ConIvsN <- na.omit(res05Ordered_ConIvsN)
upreg_ConIvsN_p0.05 <- data_ConIvsN[(as.matrix(data_ConIvsN[,6]) < 0.05)&(as.matrix(data_ConIvsN[,2]) > 0),]
dim(upreg_ConIvsN_p0.05)

## intersecting upreg genes in ConI and upreg genes in ConJ 

shared_upreg_genes_IandJ <- intersect(rownames(upreg_ConJvs18Z_p0.05),rownames(upreg_ConIvsN_p0.05))
length(shared_upreg_genes_IandJ)

library(VennDiagram)
grid.newpage()
draw.pairwise.venn(1166, 1251, 510, category = c("upreg_ConJvs18Z_p0.05", "upreg_ConIvsN_p0.05"), fill = c("light blue", "orange"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))

pdf(".../Venn_diag_intersectgns_upreg_ConI_ConJ.pdf")
library(VennDiagram)
grid.newpage()
draw.pairwise.venn(1166, 1251, 510, category = c("upreg_ConJvs18Z_p0.05", "upreg_ConIvsN_p0.05"), fill = c("light blue", "orange"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
dev.off()

write.csv(as.data.frame(shared_upreg_genes_IandJ), file=".../shared_upreg_genes_p0.05_IandJ_for_ConJvs18Z_ConIvsN.csv")


### Analysis for figure involving 3D culture of iN cells or their co-culture with mouse-astro ConC-D-I-J

dds_sub_ConCDIJ <- dds[, dds$condition %in% c("ConC", "ConD", "ConI", "ConJ")]

dds_sub_ConCDIJ$condition <- droplevels(dds_sub_ConCDIJ$condition)

rld_sub_ConCDIJ <- rlog(dds_sub_ConCDIJ, blind=FALSE)

plotPCA(rld_sub_ConCDIJ)

(data_ConCDIJ <- plotPCA(rld_sub_ConCDIJ, returnData=TRUE))

percentVar <- round(100 * attr(data_ConCDIJ, "percentVar"))

library("ggplot2")

shapes <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)

ggplot(data_ConCDIJ, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())

pdf(".../Con_C_D_I_J/PCA_Con_C_D_I_J.pdf")
ggplot(data_ConCDIJ, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
dev.off()

## differential expression

dds_ConJvsI <- dds[, dds$condition %in% c("ConI", "ConJ")]
dds_ConJvsI$condition <- relevel(dds_ConJvsI$condition, ref="ConI")
dds_ConJvsI$condition <- droplevels(dds_ConJvsI$condition)
dds_ConJvsI <- DESeq(dds_ConJvsI)

res05_ConJvsI <- results(dds_ConJvsI, alpha=0.05) 
res05Ordered_ConJvsI <- res05_ConJvsI[order(res05_ConJvsI$padj),]

write.csv(as.data.frame(res05Ordered_ConJvsI), file=".../ConJvsConI_diffexp_adjp005.csv")

## upregulated genes in ConJ for ConJvsI (with adjp<0.05)

data_ConJvsI <- na.omit(res05Ordered_ConJvsI)
upreg_ConJvsI_p0.05 <- data_ConJvsI[(as.matrix(data_ConJvsI[,6]) < 0.05)&(as.matrix(data_ConJvsI[,2]) > 0),]
dim(upreg_ConJvsI_p0.05)

### neuronal developmental and axon guidance pathways

forebrain=read.table(file=".../GO_FOREBRAIN_DEVELOPMENT.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
forebrain_J <- intersect(rownames(upreg_ConJvsI_p0.05),rownames(forebrain))
length(forebrain_J)

neuron_dev=read.table(file=".../GO_NEURON_DEVELOPMENT.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
neuron_dev_J <- intersect(rownames(upreg_ConJvsI_p0.05),rownames(neuron_dev))
length(neuron_dev_J)

cent_nervsys_dev=read.table(file=".../GO_CENTRAL_NERVOUS_SYSTEM_DEVELOPMENT.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
cent_nervsys_dev_J <- intersect(rownames(upreg_ConJvsI_p0.05),rownames(cent_nervsys_dev))
length(cent_nervsys_dev_J)

neurogenesis=read.table(file=".../GO_0022008_neurogenesis.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
neurogenesis_J <- intersect(rownames(upreg_ConJvsI_p0.05),rownames(neurogenesis))
length(neurogenesis_J)

### heatmaps for intersected genes

dds_sub_I_J <- dds[, dds$condition %in% c("ConI", "ConJ")]
dds_sub_I_J$condition <- droplevels(dds_sub_I_J$condition)

rld_sub_I_J <- rlog(dds_sub_I_J, blind=FALSE)

library("pheatmap")
library("RColorBrewer")
df <- as.data.frame(colData(dds)[c("condition")])

mt_forebrain <- assay(rld_sub_I_J)[forebrain_J,]
mt_forebrain <- mt_forebrain - rowMeans(mt_forebrain)
pheatmap(mt_forebrain, annotation_col=df, cluster_cols=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))

pdf(".../forebrain_genes_upregJvsI_adjp0.05.pdf")
pheatmap(mt_forebrain, annotation_col=df, cluster_cols=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))
dev.off()

mt_neurondev <- assay(rld_sub_I_J)[neuron_dev_J,]
mt_neurondev <- mt_neurondev - rowMeans(mt_neurondev)
pheatmap(mt_neurondev, annotation_col=df, cluster_cols=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))

pdf(".../neurondev_genes_upregJvsI_adjp0.05.pdf")
pheatmap(mt_neurondev, annotation_col=df, cluster_cols=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))
dev.off()

mt_cent_nervsys_dev <- assay(rld_sub_I_J)[cent_nervsys_dev_J,]
mt_cent_nervsys_dev <- mt_cent_nervsys_dev - rowMeans(mt_cent_nervsys_dev)
pheatmap(mt_cent_nervsys_dev, annotation_col=df, cluster_cols=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))

pdf(".../cent_nervsys_dev_genes_upregJvsI_adjp0.05.pdf")
pheatmap(mt_cent_nervsys_dev, annotation_col=df, cluster_cols=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))
dev.off()

mt_neurogenesis <- assay(rld_sub_I_J)[neurogenesis_J,]
mt_neurogenesis <- mt_neurogenesis - rowMeans(mt_neurogenesis)
pheatmap(mt_neurogenesis, annotation_col=df, cluster_cols=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))

pdf(".../neurogenesis_genes_upregJvsI_adjp0.05.pdf")
pheatmap(mt_neurogenesis, annotation_col=df, cluster_cols=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))
dev.off()
