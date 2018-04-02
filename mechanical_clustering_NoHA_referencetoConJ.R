
countData = read.table(file=".../rsem.genes.counts.matrix",header=TRUE,sep="\t",row.names=1,check.names = FALSE)
colData = read.table(file=".../samples_described.txt",head=TRUE,sep="\t",row.names=1,check.names = FALSE)
rownames(colData)=sub("-","_",rownames(colData), fixed=T)
countData <- countData[ , rownames(colData)]
countData = round(countData)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)


dds_3.125Ca <- dds[, dds$condition %in% c("ConJ", "Con14X")]
dds_3.125Ca$condition <- relevel(dds_3.125Ca$condition, ref="ConJ")
dds_3.125Ca$condition <- droplevels(dds_3.125Ca$condition)
dds_3.125Ca <- DESeq(dds_3.125Ca)

res05_3.125Ca <- results(dds_3.125Ca, alpha=0.05)
res05Ordered_3.125Ca <- res05_3.125Ca[order(res05_3.125Ca$padj),]
data_3.125Ca <- na.omit(res05Ordered_3.125Ca)
dim(data_3.125Ca)
fltr_3.125Ca <- data_3.125Ca[(as.matrix(data_3.125Ca[,6]) < 0.01)&(as.matrix(data_3.125Ca[,2]) > 1 | as.matrix(data_3.125Ca[,2]) < -1),]
dim(fltr_3.125Ca)

dds_6.25Ca <- dds[, dds$condition %in% c("ConJ", "Con4X")]
dds_6.25Ca$condition <- droplevels(dds_6.25Ca$condition)
dds_6.25Ca$condition <- relevel(dds_6.25Ca$condition, ref="ConJ")
dds_6.25Ca$condition <- droplevels(dds_6.25Ca$condition)
dds_6.25Ca <- DESeq(dds_6.25Ca)

res05_6.25Ca <- results(dds_6.25Ca, alpha=0.05)
res05Ordered_6.25Ca <- res05_6.25Ca[order(res05_6.25Ca$padj),]
data_6.25Ca <- na.omit(res05Ordered_6.25Ca)
dim(data_6.25Ca)
fltr_6.25Ca <- data_6.25Ca[(as.matrix(data_6.25Ca[,6]) < 0.01)&(as.matrix(data_6.25Ca[,2]) > 1 | as.matrix(data_6.25Ca[,2]) < -1),] 
dim(fltr_6.25Ca)

dds_12.5Ca <- dds[, dds$condition %in% c("ConJ", "Con1X")]
dds_12.5Ca$condition <- relevel(dds_12.5Ca$condition, ref="ConJ")
dds_12.5Ca$condition <- droplevels(dds_12.5Ca$condition)
dds_12.5Ca <- DESeq(dds_12.5Ca)

res05_12.5Ca <- results(dds_12.5Ca, alpha=0.05)
res05Ordered_12.5Ca <- res05_12.5Ca[order(res05_12.5Ca$padj),]
data_12.5Ca <- na.omit(res05Ordered_12.5Ca)
dim(data_12.5Ca)
fltr_12.5Ca <- data_12.5Ca[(as.matrix(data_12.5Ca[,6]) < 0.01)&(as.matrix(data_12.5Ca[,2]) > 1 | as.matrix(data_12.5Ca[,2]) < -1),]
dim(fltr_12.5Ca)

dds_25Ca <- dds[, dds$condition %in% c("ConJ", "ConT")]
dds_25Ca$condition <- relevel(dds_25Ca$condition, ref="ConJ")
dds_25Ca$condition <- droplevels(dds_25Ca$condition)
dds_25Ca <- DESeq(dds_25Ca)

res05_25Ca <- results(dds_25Ca, alpha=0.05)
res05Ordered_25Ca <- res05_25Ca[order(res05_25Ca$padj),]
data_25Ca <- na.omit(res05Ordered_25Ca)
dim(data_25Ca)
fltr_25Ca <- data_25Ca[(as.matrix(data_25Ca[,6]) < 0.01)&(as.matrix(data_25Ca[,2]) > 1 | as.matrix(data_25Ca[,2]) < -1),]
dim(fltr_25Ca)

rows_3.125Ca <- rownames(fltr_3.125Ca)
rows_6.25Ca <- rownames(fltr_6.25Ca)
combine_3.125Ca6.25Ca <- c(rows_3.125Ca, rows_6.25Ca)
genes_3.125Ca6.25Ca <- unique(combine_3.125Ca6.25Ca)
length(genes_3.125Ca6.25Ca)

rows_12.5Ca <- rownames(fltr_12.5Ca)
rows_25Ca <- rownames(fltr_25Ca)
combine_12.5Ca25Ca <- c(rows_12.5Ca, rows_25Ca)
genes_12.5Ca25Ca <- unique(combine_12.5Ca25Ca)
length(genes_12.5Ca25Ca)

combine_all <- c(genes_3.125Ca6.25Ca, genes_12.5Ca25Ca)
length(combine_all)
mech_genes <- unique(combine_all)
length(mech_genes)

dds_sub <- dds[, dds$condition %in% c("ConJ", "Con14X", "Con4X", "Con1X", "ConT")]
dds_sub$condition <- droplevels(dds_sub$condition)
rld_sub <- rlog(dds_sub, blind=FALSE)

plotPCA(rld_sub)
(data <- plotPCA(rld_sub, returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))
library("ggplot2")
#shapes <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)
#shapes <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,25)
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())

pdf(".../PCA_ConJ_14X_4X_1X_T.pdf")
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
dev.off()

library("pheatmap")
library("RColorBrewer")
df <- as.data.frame(colData(dds_sub)[c("condition")])

mat <- assay(rld_sub)[mech_genes,]
mat <- mat - rowMeans(mat)

km<- kmeans(mat,4)
m.kmeans<- cbind(mat, km$cluster)
dim(m.kmeans)
o<- order(m.kmeans[,16])
m.kmeans<- m.kmeans[o,]

a<-m.kmeans[,1:15]
a <- a[,c("10A_S28","10B_S29","10C_S30","14A060216_S40", "14B060216_S41", "14C060216_S42", "4A060216_S10" , "4B060216_S11",  "4C060216_S12","1A060116_S1" ,  "1B060116_S2" , "1C060116_S3",  "19A_S55","19B_S56","19C_S57")]
pheatmap( a, cluster_rows = F, cluster_cols = F,  annotation_col=df, show_rownames=FALSE, color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))

pdf(".../kmeans_4clusters.pdf")
pheatmap( a, cluster_rows = F, cluster_cols = F,  annotation_col=df, show_rownames=FALSE, color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()

write.csv(as.data.frame(a), file=".../kmeans_4clusters_columnsorted.csv")
write.csv(as.data.frame(m.kmeans), file=".../kmeans_4clusters_withclusternumbers.csv")



forebrain=read.table(file=".../GO_FOREBRAIN_DEVELOPMENT.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
forebrain_mech <- intersect(rownames(a),rownames(forebrain))
length(forebrain_mech)
mt <- assay(rld_sub)[forebrain_mech,]
mt <- mt - rowMeans(mt)

mt <- mt[,c("10A_S28","10B_S29","10C_S30","14A060216_S40", "14B060216_S41", "14C060216_S42", "4A060216_S10" , "4B060216_S11",  "4C060216_S12","1A060116_S1" ,  "1B060116_S2" , "1C060116_S3",  "19A_S55","19B_S56","19C_S57")]
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6.5,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))

pdf(".../forebrain_dev_mechgenes_intersect.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6.5,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()

kegg_axon_guid=read.table(file=".../kegg_axon_guidance_geneset.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
axonguid_mech <- intersect(rownames(a),rownames(kegg_axon_guid))
length(axonguid_mech)
mt <- assay(rld_sub)[axonguid_mech,]
mt <- mt - rowMeans(mt)

mt <- mt[,c("10A_S28","10B_S29","10C_S30","14A060216_S40", "14B060216_S41", "14C060216_S42", "4A060216_S10" , "4B060216_S11",  "4C060216_S12","1A060116_S1" ,  "1B060116_S2" , "1C060116_S3",  "19A_S55","19B_S56","19C_S57")]
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6.5,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))

pdf(".../axon_guidance_kegg_mechgenes_intersect.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6.5,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()

cent_nervsys_dev=read.table(file=".../GO_CENTRAL_NERVOUS_SYSTEM_DEVELOPMENT.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
cent_nrv_sys_mech <- intersect(rownames(a),rownames(cent_nervsys_dev))
length(cent_nrv_sys_mech)
mt <- assay(rld_sub)[cent_nrv_sys_mech,]
mt <- mt - rowMeans(mt)

mt <- mt[,c("10A_S28","10B_S29","10C_S30","14A060216_S40", "14B060216_S41", "14C060216_S42", "4A060216_S10" , "4B060216_S11",  "4C060216_S12","1A060116_S1" ,  "1B060116_S2" , "1C060116_S3",  "19A_S55","19B_S56","19C_S57")]
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=3.5,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))

pdf(".../cent_nerv_sys_dev_mechgenes_intersect.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=3.5,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()

neuron_dev=read.table(file=".../GO_NEURON_DEVELOPMENT.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
neurondev_mech <- intersect(rownames(a),rownames(neuron_dev))
length(neurondev_mech)
mt <- assay(rld_sub)[neurondev_mech,]
mt <- mt - rowMeans(mt)

mt <- mt[,c("10A_S28","10B_S29","10C_S30","14A060216_S40", "14B060216_S41", "14C060216_S42", "4A060216_S10" , "4B060216_S11",  "4C060216_S12","1A060116_S1" ,  "1B060116_S2" , "1C060116_S3",  "19A_S55","19B_S56","19C_S57")]
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))

pdf(".../neuron_dev_mechgenes_intersect.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()


autism_syndromic_genes <- read.table(file=".../Autism_genes.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
autsm_synd_mech_gns <- intersect(rownames(a),rownames(autism_syndromic_genes))
length(autsm_synd_mech_gns)
mt <- assay(rld_sub)[autsm_synd_mech_gns,]
mt <- mt - rowMeans(mt)

mt <- mt[,c("10A_S28","10B_S29","10C_S30","14A060216_S40", "14B060216_S41", "14C060216_S42", "4A060216_S10" , "4B060216_S11",  "4C060216_S12","1A060116_S1" ,  "1B060116_S2" , "1C060116_S3",  "19A_S55","19B_S56","19C_S57")]
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))

pdf(".../autsm_synd_mech_gns_intersect.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()


autism_highconf_genes <- read.table(file=".../autism_highconfidence_genes.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
autsm_highconf_mech_gns <- intersect(rownames(a),rownames(autism_highconf_genes))
length(autsm_highconf_mech_gns)
mt <- assay(rld_sub)[autsm_highconf_mech_gns,]
mt <- mt - rowMeans(mt)

mt <- mt[,c("10A_S28","10B_S29","10C_S30","14A060216_S40", "14B060216_S41", "14C060216_S42", "4A060216_S10" , "4B060216_S11",  "4C060216_S12","1A060116_S1" ,  "1B060116_S2" , "1C060116_S3",  "19A_S55","19B_S56","19C_S57")]
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))

pdf(".../autsm_highconf_mech_gns_intersect.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()


Alsod_gene_Data <- read.table(file=".../Alsod_gene_Data.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
Alsod_gene_Data_mech_gns <- intersect(rownames(a),rownames(Alsod_gene_Data))
length(Alsod_gene_Data_mech_gns)
mt <- assay(rld_sub)[Alsod_gene_Data_mech_gns,]
mt <- mt - rowMeans(mt)

mt <- mt[,c("10A_S28","10B_S29","10C_S30","14A060216_S40", "14B060216_S41", "14C060216_S42", "4A060216_S10" , "4B060216_S11",  "4C060216_S12","1A060116_S1" ,  "1B060116_S2" , "1C060116_S3",  "19A_S55","19B_S56","19C_S57")]
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))

pdf(".../Alsod_gene_Data_mech_gns_intersect.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()


ALSgenesorg_topgenes <- read.table(file=".../ALSgenesorg_topgenes.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
ALSgenesorg_topgenes_mech_gns <- intersect(rownames(a),rownames(ALSgenesorg_topgenes))
length(ALSgenesorg_topgenes_mech_gns)
mt <- assay(rld_sub)[ALSgenesorg_topgenes_mech_gns,]
mt <- mt - rowMeans(mt)

mt <- mt[,c("10A_S28","10B_S29","10C_S30","14A060216_S40", "14B060216_S41", "14C060216_S42", "4A060216_S10" , "4B060216_S11",  "4C060216_S12","1A060116_S1" ,  "1B060116_S2" , "1C060116_S3",  "19A_S55","19B_S56","19C_S57")]
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))

pdf(".../ALSgenesorg_topgenes_mech_gns_intersect.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()


Alzheimer_genes <- read.table(file=".../Alzheimer_genes.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
Alzheimer_genes_mech_gns <- intersect(rownames(a),rownames(Alzheimer_genes))
length(Alzheimer_genes_mech_gns)
mt <- assay(rld_sub)[Alzheimer_genes_mech_gns,]
mt <- mt - rowMeans(mt)

mt <- mt[,c("10A_S28","10B_S29","10C_S30","14A060216_S40", "14B060216_S41", "14C060216_S42", "4A060216_S10" , "4B060216_S11",  "4C060216_S12","1A060116_S1" ,  "1B060116_S2" , "1C060116_S3",  "19A_S55","19B_S56","19C_S57")]
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))

pdf(".../Alzheimer_genes_mech_gns_intersect.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()


Parkinson_genes <- read.table(file=".../Parkinson_genes.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
Parkinson_genes_mech_gns <- intersect(rownames(a),rownames(Parkinson_genes))
length(Parkinson_genes_mech_gns)
mt <- assay(rld_sub)[Parkinson_genes_mech_gns,]
mt <- mt - rowMeans(mt)

mt <- mt[,c("10A_S28","10B_S29","10C_S30","14A060216_S40", "14B060216_S41", "14C060216_S42", "4A060216_S10" , "4B060216_S11",  "4C060216_S12","1A060116_S1" ,  "1B060116_S2" , "1C060116_S3",  "19A_S55","19B_S56","19C_S57")]
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))

pdf(".../Parkinson_genes_mech_gns_intersect.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()



### heatmap of selected neuronal genes for main figure 

selected <- c("DLG4","PAK2","ATF5","CLU","NDST1","NEUROD1","RELN","NFIB","EPHA3","SEMA3C","SOD1","ETV1","PTN","RAC3","EPHB1","NF1","LHX9","SLC1A2","ARHGEF12","SECISBP2")
mt <- assay(rld_sub)[selected,]
mt <- mt - rowMeans(mt)

mt <- mt[,c("10A_S28","10B_S29","10C_S30","14A060216_S40", "14B060216_S41", "14C060216_S42", "4A060216_S10" , "4B060216_S11",  "4C060216_S12","1A060116_S1" ,  "1B060116_S2" , "1C060116_S3",  "19A_S55","19B_S56","19C_S57")]
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=7,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))

pdf(".../selected_neuronal_genes_for_main_figure.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=7,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()



### combined heatmap of selected neuronal genes and disease associated genes for main figure 

selected <- c("DLG4","PAK2","ATF5","CLU","NDST1","NEUROD1","RELN","NFIB","EPHA3","SEMA3C","SOD1","ETV1","PTN","RAC3","EPHB1","NF1","LHX9","SLC1A2","ARHGEF12","SECISBP2","DHCR7","ADSL","BCL6","AGT","CACNA1C","RELN","SEMA6A","TRPM7","C12orf57","SOD1","CST3","ZSWIM7","ASXL3","MBD5","SLC1A2","CRIM1")
mt <- assay(rld_sub)[selected,]
mt <- mt - rowMeans(mt)

mt <- mt[,c("10A_S28","10B_S29","10C_S30","14A060216_S40", "14B060216_S41", "14C060216_S42", "4A060216_S10" , "4B060216_S11",  "4C060216_S12","1A060116_S1" ,  "1B060116_S2" , "1C060116_S3",  "19A_S55","19B_S56","19C_S57")]
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=7,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))

pdf(".../selected_combined_neuronal_genes_&_disease_assctd_genes_for_main_figure.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=7,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()



