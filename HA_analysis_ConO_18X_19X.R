
countData = read.table(file=".../rsem.genes.counts.matrix",header=TRUE,sep="\t",row.names=1,check.names = FALSE)
colData = read.table(file=".../samples_described.txt",head=TRUE,sep="\t",row.names=1,check.names = FALSE)

rownames(colData)=sub("-","_",rownames(colData), fixed=T)
countData <- countData[ , rownames(colData)]
countData = round(countData)
 
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

dds_sub <- dds[, dds$condition %in% c("ConO", "Con18X", "Con19X")]
dds_sub$condition <- droplevels(dds_sub$condition)

rld_sub <- rlog(dds_sub, blind=FALSE)

### plotting PCA ####

plotPCA(rld_sub)
(data <- plotPCA(rld_sub, returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))
library("ggplot2")
#shapes <- c(15,7,17)
shapes <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,25)

ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())

pdf(".../PCA_ConO_18X_19X.pdf")
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4.5) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=5) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=6) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
dev.off()

### Differential Expression Analysis ###

dds_1HA <- dds[, dds$condition %in% c("ConO", "Con18X")]
dds_1HA$condition <- relevel(dds_1HA$condition, ref="ConO")
dds_1HA$condition <- droplevels(dds_1HA$condition)
dds_1HA <- DESeq(dds_1HA)

res05_1HA <- results(dds_1HA, alpha=0.05) 
res05Ordered_1HA <- res05_1HA[order(res05_1HA$padj),]

write.csv(as.data.frame(res05Ordered_1HA), file=".../Con18XvsConO_diffexp_adjp005.csv")

dds_1.5HA <- dds[, dds$condition %in% c("ConO", "Con19X")]
dds_1.5HA$condition <- relevel(dds_1.5HA$condition, ref="ConO")
dds_1.5HA$condition <- droplevels(dds_1.5HA$condition)
dds_1.5HA <- DESeq(dds_1.5HA)

res05_1.5HA <- results(dds_1.5HA, alpha=0.05)
res05Ordered_1.5HA <- res05_1.5HA[order(res05_1.5HA$padj),]

write.csv(as.data.frame(res05Ordered_1.5HA), file=".../Con19XvsConO_diffexp_adjp005.csv")

### Clustering for diff exp genes with p<0.01 and log2foldchange >1 or <-1 ####

data_1HA <- na.omit(res05Ordered_1HA)
fltr_1HA <- data_1HA[(as.matrix(data_1HA[,6]) < 0.01)&(as.matrix(data_1HA[,2]) > 1 | as.matrix(data_1HA[,2]) < -1),]
dim(fltr_1HA)

data_1.5HA <- na.omit(res05Ordered_1.5HA)
fltr_1.5HA <- data_1.5HA[(as.matrix(data_1.5HA[,6]) < 0.01)&(as.matrix(data_1.5HA[,2]) > 1 | as.matrix(data_1.5HA[,2]) < -1),]
dim(fltr_1.5HA)

rows_1HA <- rownames(fltr_1HA)
rows_1.5HA <- rownames(fltr_1.5HA)
combine_HA <- c(rows_1HA,rows_1.5HA)

HA_genes <- unique(combine_HA)
length(HA_genes)

mat <- assay(rld_sub)[HA_genes,]
mat <- mat - rowMeans(mat)

library("pheatmap")
library("RColorBrewer")
df <- as.data.frame(colData(dds)[c("condition")])

km<- kmeans(mat,2)
m.kmeans<- cbind(mat, km$cluster)
dim(m.kmeans)

o<- order(m.kmeans[,10])
m.kmeans<- m.kmeans[o,]

pheatmap( m.kmeans[,1:9],cluster_rows = F, cluster_cols = F,  annotation_col=df, show_rownames=FALSE, color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))

pdf(".../HA_kmeans_2clusters_ConO_18X_19X.pdf")
pheatmap( m.kmeans[,1:9],cluster_rows = F, cluster_cols = F,  annotation_col=df, show_rownames=FALSE, color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()

write.csv(as.data.frame(m.kmeans), file=".../ConO_Con18X_Con19X_kmeans2clusters.csv")


#### Intersecting Upregulated genes in Con18XvsConO (with adjp<0.05) with gene sets ###

upreg_18XvsO_p0.05 <- data_1HA[(as.matrix(data_1HA[,6]) < 0.05)&(as.matrix(data_1HA[,2]) > 0),]

### neuronal developmental and axon guidance pathways

forebrain=read.table(file=".../GO_FOREBRAIN_DEVELOPMENT.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
forebrain_18X <- intersect(rownames(upreg_18XvsO_p0.05),rownames(forebrain))

neuron_dev=read.table(file=".../GO_NEURON_DEVELOPMENT.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
neuron_dev_18X <- intersect(rownames(upreg_18XvsO_p0.05),rownames(neuron_dev))

kegg_axon_guid=read.table(file=".../kegg_axon_guidance_geneset.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
axon_guid_18X <- intersect(rownames(upreg_18XvsO_p0.05),rownames(kegg_axon_guid))

cent_nervsys_dev=read.table(file=".../GO_CENTRAL_NERVOUS_SYSTEM_DEVELOPMENT.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
cent_nervsys_dev_18X <- intersect(rownames(upreg_18XvsO_p0.05),rownames(cent_nervsys_dev))

reactome_axon_guid=read.table(file=".../REACTOME_AXON_GUIDANCE_REACT_18266.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
reactome_axon_guid_18X <- intersect(rownames(upreg_18XvsO_p0.05),rownames(reactome_axon_guid))
length(reactome_axon_guid_18X)


### channel activity 

volt_gated_ion_ch=read.table(file=".../GO_VOLTAGE_GATED_ION_CHNNL_ACT.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
volt_gated_ion_ch_18X <- intersect(rownames(upreg_18XvsO_p0.05),rownames(volt_gated_ion_ch))
length(volt_gated_ion_ch_18X)

ion_gated_ch_act=read.table(file=".../GO_ION_GATED_CHANNEL_ACTIVITY.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
ion_gated_ch_act_18X <- intersect(rownames(upreg_18XvsO_p0.05),rownames(ion_gated_ch_act))
length(ion_gated_ch_act_18X)

cal_ch_cmplx=read.table(file=".../GO_CALCIUM_CHANNEL_COMPLEX.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
cal_ch_cmplx_18X <- intersect(rownames(upreg_18XvsO_p0.05),rownames(cal_ch_cmplx))
length(cal_ch_cmplx_18X)

chlrd_ch_cmplx=read.table(file=".../GO_CHLORIDE_CHANNEL_COMPLEX.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
chlrd_ch_cmplx_18X <- intersect(rownames(upreg_18XvsO_p0.05),rownames(chlrd_ch_cmplx))
length(chlrd_ch_cmplx_18X)

sodium_ch_cmplx=read.table(file=".../sodium_channel_complex.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
sodium_ch_cmplx_18X <- intersect(rownames(upreg_18XvsO_p0.05),rownames(sodium_ch_cmplx))
length(sodium_ch_cmplx_18X)

sodium_ch_act=read.table(file=".../GO_SODIUM_CHANNEL_ACTIVITY.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
sodium_ch_act_18X <- intersect(rownames(upreg_18XvsO_p0.05),rownames(sodium_ch_act))
length(sodium_ch_act_18X)

potsm_ch_act=read.table(file=".../Potassium_Chnl_Act.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
potsm_ch_act_18X <- intersect(rownames(upreg_18XvsO_p0.05),rownames(potsm_ch_act))
length(potsm_ch_act_18X)


### heatmaps for intersected genes

dds_sub_18X_O <- dds[, dds$condition %in% c("ConO", "Con18X")]
dds_sub_18X_O$condition <- droplevels(dds_sub_18X_O$condition)

rld_sub_18X_O <- rlog(dds_sub_18X_O, blind=FALSE)

mt_forebrain <- assay(rld_sub_18X_O)[forebrain_18X,]
mt_forebrain <- mt_forebrain - rowMeans(mt_forebrain)
pheatmap(mt_forebrain, annotation_col=df, cluster_cols=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))

pdf(".../forebrain_genes_upreg18XvsO_adjp0.05.pdf")
pheatmap(mt_forebrain, annotation_col=df, cluster_cols=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))
dev.off()

mt_neurondev <- assay(rld_sub_18X_O)[neuron_dev_18X,]
mt_neurondev <- mt_neurondev - rowMeans(mt_neurondev)
pheatmap(mt_neurondev, annotation_col=df, cluster_cols=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))

pdf(".../neurondev_genes_upreg18XvsO_adjp0.05.pdf")
pheatmap(mt_neurondev, annotation_col=df, cluster_cols=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))
dev.off()

mt_axon_guid <- assay(rld_sub_18X_O)[axon_guid_18X,]
mt_axon_guid <- mt_axon_guid - rowMeans(mt_axon_guid)
pheatmap(mt_axon_guid, annotation_col=df, cluster_cols=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))

pdf(".../kegg_axon_guid_genes_upreg18XvsO_adjp0.05.pdf")
pheatmap(mt_axon_guid, annotation_col=df, cluster_cols=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))
dev.off()

mt_cent_nervsys_dev <- assay(rld_sub_18X_O)[cent_nervsys_dev_18X,]
mt_cent_nervsys_dev <- mt_cent_nervsys_dev - rowMeans(mt_cent_nervsys_dev)
pheatmap(mt_cent_nervsys_dev, annotation_col=df, cluster_cols=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))

pdf(".../cent_nervsys_dev_genes_upreg18XvsO_adjp0.05.pdf")
pheatmap(mt_cent_nervsys_dev, annotation_col=df, cluster_cols=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))
dev.off()

mt_reactome_axon_guid <- assay(rld_sub_18X_O)[reactome_axon_guid_18X,]
mt_reactome_axon_guid <- mt_reactome_axon_guid - rowMeans(mt_reactome_axon_guid)
pheatmap(mt_reactome_axon_guid, annotation_col=df, cluster_cols=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))

pdf(".../reactome_axon_guid_genes_upreg18XvsO_adjp0.05.pdf")
pheatmap(mt_reactome_axon_guid, annotation_col=df, cluster_cols=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))
dev.off()


### heatmap for selected genes from intersected genes

selected <- c( "NUMB", "DLG4", "NEUROD4", "UBB", "CLU", "BAG3", "NEFH", "EGFR", "IGF2BP1","ATM", "EPOR", "ARF4", "NHLH1", "KCNN1", "ABCC9")

neuronal_genes_18X <- assay(rld_sub_18X_O)[selected,]
neuronal_genes_18X <- neuronal_genes_18X - rowMeans(neuronal_genes_18X)
pheatmap(neuronal_genes_18X, annotation_col=df, cluster_cols=FALSE, cluster_rows=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))

pdf(".../selected_genes_mainfigure_upreg18XvsO_adjp0.05.pdf")
pheatmap(neuronal_genes_18X, annotation_col=df, cluster_cols=FALSE, cluster_rows=FALSE, fontsize=10,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))
dev.off()
              
