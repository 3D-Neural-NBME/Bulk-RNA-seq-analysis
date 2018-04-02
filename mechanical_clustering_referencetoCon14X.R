
countData = read.table(file=".../rsem.genes.counts.matrix",header=TRUE,sep="\t",row.names=1,check.names = FALSE)
colData = read.table(file=".../samples_described.txt",head=TRUE,sep="\t",row.names=1,check.names = FALSE)
rownames(colData)=sub("-","_",rownames(colData), fixed=T)
countData <- countData[ , rownames(colData)]
countData = round(countData)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)


####### using DE genes with for p<0.05 and log2(fold change) > 0.75 or log2(fold change) < -0.75 #############

dds_6.25Ca <- dds[, dds$condition %in% c("Con14X", "Con4X")]
dds_6.25Ca$condition <- relevel(dds_6.25Ca$condition, ref="Con14X")
dds_6.25Ca$condition <- droplevels(dds_6.25Ca$condition)
dds_6.25Ca <- DESeq(dds_6.25Ca)
res05_6.25Ca <- results(dds_6.25Ca, alpha=0.05)
res05Ordered_6.25Ca <- res05_6.25Ca[order(res05_6.25Ca$padj),]
data_6.25Ca <- na.omit(res05Ordered_6.25Ca)
dim(data_6.25Ca)
fltr_6.25Ca <- data_6.25Ca[(as.matrix(data_6.25Ca[,6]) < 0.05)&(as.matrix(data_6.25Ca[,2]) > 0.75 | as.matrix(data_6.25Ca[,2]) < -0.75),]
dim(fltr_6.25Ca)

dds_12.5Ca <- dds[, dds$condition %in% c("Con14X", "Con1X")]
dds_12.5Ca$condition <- relevel(dds_12.5Ca$condition, ref="Con14X")
dds_12.5Ca$condition <- droplevels(dds_12.5Ca$condition)
dds_12.5Ca <- DESeq(dds_12.5Ca)
res05_12.5Ca <- results(dds_12.5Ca, alpha=0.05)
res05Ordered_12.5Ca <- res05_12.5Ca[order(res05_12.5Ca$padj),]
data_12.5Ca <- na.omit(res05Ordered_12.5Ca)
dim(data_12.5Ca)
fltr_12.5Ca <- data_12.5Ca[(as.matrix(data_12.5Ca[,6]) < 0.05)&(as.matrix(data_12.5Ca[,2]) > 0.75 | as.matrix(data_12.5Ca[,2]) < -0.75),]
dim(fltr_12.5Ca)


dds_25Ca <- dds[, dds$condition %in% c("Con14X", "ConT")]
dds_25Ca$condition <- relevel(dds_25Ca$condition, ref="Con14X")
dds_25Ca$condition <- droplevels(dds_25Ca$condition)
dds_25Ca <- DESeq(dds_25Ca)
res05_25Ca <- results(dds_25Ca, alpha=0.05)
res05Ordered_25Ca <- res05_25Ca[order(res05_25Ca$padj),]
data_25Ca <- na.omit(res05Ordered_25Ca)
dim(data_25Ca)
fltr_25Ca <- data_25Ca[(as.matrix(data_25Ca[,6]) < 0.05)&(as.matrix(data_25Ca[,2]) > 0.75 | as.matrix(data_25Ca[,2]) < -0.75),]
dim(fltr_25Ca)

rows_6.25Ca <- rownames(fltr_6.25Ca)
rows_12.5Ca <- rownames(fltr_12.5Ca)
combine_6.25Ca12.5Ca <- c(rows_6.25Ca, rows_12.5Ca)
genes_6.25Ca12.5Ca <- unique(combine_6.25Ca12.5Ca)
length(genes_6.25Ca12.5Ca)
rows_25Ca <- rownames(fltr_25Ca)
combine_all <- c(combine_6.25Ca12.5Ca,rows_25Ca)
genes_all <- unique(combine_all)
length(genes_all)
mech_genes <- genes_all

dds_sub <- dds[, dds$condition %in% c("Con14X", "Con4X", "Con1X", "ConT")]
dds_sub$condition <- droplevels(dds_sub$condition)
rld_sub <- rlog(dds_sub, blind=FALSE)

library("pheatmap")
library("RColorBrewer")
df <- as.data.frame(colData(dds_sub)[c("condition")])

mat <- assay(rld_sub)[mech_genes,]
mat <- mat - rowMeans(mat)



km<- kmeans(mat,4)
m.kmeans<- cbind(mat, km$cluster)
dim(m.kmeans)
o<- order(m.kmeans[,13])
m.kmeans<- m.kmeans[o,]
a<-m.kmeans[,1:12]
a <- a[,c("14A060216_S40", "14B060216_S41", "14C060216_S42", "4A060216_S10" , "4B060216_S11",  "4C060216_S12","1A060116_S1" ,  "1B060116_S2" , "1C060116_S3",  "19A_S55","19B_S56","19C_S57")]
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
mt <- mt[,c("14A060216_S40", "14B060216_S41", "14C060216_S42", "4A060216_S10" , "4B060216_S11",  "4C060216_S12","1A060116_S1" ,  "1B060116_S2" , "1C060116_S3",  "19A_S55","19B_S56","19C_S57")]
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6.5,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
pdf(".../forebrain_dev_mechgenes_intersect.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6.5,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()


kegg_axon_guid=read.table(file=".../kegg_axon_guidance_geneset.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
axonguid_mech <- intersect(rownames(a),rownames(kegg_axon_guid))
length(axonguid_mech)
mt <- assay(rld_sub)[axonguid_mech,]
mt <- mt - rowMeans(mt)
mt <- mt[,c("14A060216_S40", "14B060216_S41", "14C060216_S42", "4A060216_S10" , "4B060216_S11",  "4C060216_S12","1A060116_S1" ,  "1B060116_S2" , "1C060116_S3",  "19A_S55","19B_S56","19C_S57")]
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6.5,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
pdf(".../axon_guidance_kegg_mechgenes_intersect.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6.5,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()


neuron_dev=read.table(file=".../GO_NEURON_DEVELOPMENT.txt",head=FALSE,sep="\t",row.names=1,check.names = FALSE)
neurondev_mech <- intersect(rownames(a),rownames(neuron_dev))
length(neurondev_mech)
mt <- assay(rld_sub)[neurondev_mech,]
mt <- mt - rowMeans(mt)
mt <- mt[,c("14A060216_S40", "14B060216_S41", "14C060216_S42", "4A060216_S10" , "4B060216_S11",  "4C060216_S12","1A060116_S1" ,  "1B060116_S2" , "1C060116_S3",  "19A_S55","19B_S56","19C_S57")]
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
pdf(".../neuron_dev_mechgenes_intersect.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=6,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()


### heatmap of selected neuronal genes for main figure 

selected <- c("BAD","ID4","ARHGEF12", "PLXNA1","GSK3B","RAC3","NF1","GRIN3A","KIRREL3","SLC4A7","SECISBP2","STMN1","PTN","GPM6A","DLG4","NFIB","ROBO1","UNC5C","SEMA3C","EPHA3","SLIT3","KIF5A","CHD7","CLU","PCDH9","EPHB1","SEMA4F","NRAS","NLGN3","CNTN1")
mt <- assay(rld_sub)[selected,]
mt <- mt - rowMeans(mt)

mt <- mt[,c("14A060216_S40", "14B060216_S41", "14C060216_S42", "4A060216_S10" , "4B060216_S11",  "4C060216_S12","1A060116_S1" ,  "1B060116_S2" , "1C060116_S3",  "19A_S55","19B_S56","19C_S57")]
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=7,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))

pdf(".../selected_neuronal_genes_for_main_figure.pdf")
pheatmap(mt, annotation_col=df, cluster_rows =FALSE, cluster_cols=FALSE, fontsize=7,color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
dev.off()
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              