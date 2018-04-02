
countData = read.table(file=".../rsem.genes.counts.matrix",header=TRUE,sep="\t",row.names=1,check.names = FALSE)
colData = read.table(file=".../samples_described.txt",head=TRUE,sep="\t",row.names=1,check.names = FALSE)
rownames(colData)=sub("-","_",rownames(colData), fixed=T)
countData <- countData[ , rownames(colData)]
countData = round(countData)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

dds_sub <- dds[, dds$condition %in% c("ConI","ConG","ConJ","ConH","ConK","ConO","Con18X","Con17X","Con3X","Con7X","Con2X","ConP","Con12X","Con10X","Con6X","Con5X","Con8X")]
dds_sub$condition <- droplevels(dds_sub$condition)
rld_sub <- rlog(dds_sub, blind=FALSE)

### z-score heatmap of selected ASD, ALS, AD and PD genes in 3D conditions for main figure ####

selected_genes <- c("CUL3","TRIP12","UBE3A","CACNA1C","ZSWIM7","SOD1","NTNG1","CLU","SLC24A4","SNCA","STK39","DLG2")
mat <- assay(rld_sub)[selected_genes,]
df <- as.data.frame(colData(dds_sub)[c("condition")])

mat <- mat[,c("9A_S25","9B_S26","9C_S27","7A_S19","7B_S20","7C_S21","10A_S28","10B_S29","10C_S30","8A_S22","8B_S23","8C_S24","11A_S31","11B_S32","11C_S33","15A_S43","15B_S44","15C_S45",
              "18A060216_S51","18B060216_S52","18C060216_S53","17A060216_S48","17B060216_S49","17C060216_S50","3A060116_S7","3B060116_S8","3C060116_S9","7A060216_S19","7B060216_S20","7C060216_S21",
              "2A060116_S4","2B060116_S5","2C060116_S6","16A_S46","16B_S47","16C_S48","12A060216_S34","12B060216_S35","12C060216_S36","10A060216_S28","10B060216_S29","10C060216_S30",
              "6A060216_S16","6B060216_S17","6C060216_S18","5A060216_S13","5B060216_S14","5C060216_S15","8A060216_S22","8B060216_S23","8C060216_S24")]

library("gplots")
library("RColorBrewer")
heatmap.2(mat,scale="row",trace="none", dendrogram="none", Rowv = FALSE,Colv = FALSE, col = colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255))
              
pdf(".../selected_disease_genes_zscore.pdf")
heatmap.2(mat,scale="row",trace="none", dendrogram="none", Rowv = FALSE,Colv = FALSE, col = colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255))
dev.off()

### z-score heatmap of selected ASD, ALS, AD and PD genes in 3D conditions for supplementary figure ####

selected_genes <- c("CHD8","SETD5","PTEN","SCN2A","MECP2","SMARCA2","C12orf57","LHFP","C9orf72","TARDBP","TBK1","CNTF","PARK7","BIN1","PICALM","FERMT2","CELF1","ABCA7","ASH1L","BCKDK","MAPT","TMEM229B","LRRK2")
mat <- assay(rld_sub)[selected_genes,]
mat <- mat[,c("9A_S25","9B_S26","9C_S27","7A_S19","7B_S20","7C_S21","10A_S28","10B_S29","10C_S30","8A_S22","8B_S23","8C_S24","11A_S31","11B_S32","11C_S33","15A_S43","15B_S44","15C_S45",
              "18A060216_S51","18B060216_S52","18C060216_S53","17A060216_S48","17B060216_S49","17C060216_S50","3A060116_S7","3B060116_S8","3C060116_S9","7A060216_S19","7B060216_S20","7C060216_S21",
              "2A060116_S4","2B060116_S5","2C060116_S6","16A_S46","16B_S47","16C_S48","12A060216_S34","12B060216_S35","12C060216_S36","10A060216_S28","10B060216_S29","10C060216_S30",
              "6A060216_S16","6B060216_S17","6C060216_S18","5A060216_S13","5B060216_S14","5C060216_S15","8A060216_S22","8B060216_S23","8C060216_S24")]

heatmap.2(mat,scale="row",trace="none", dendrogram="none", Rowv = FALSE,Colv = FALSE, col = colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255))

pdf(".../for_supp_fig_selected_disease_genes_zscore.pdf")
heatmap.2(mat,scale="row",trace="none", dendrogram="none", Rowv = FALSE,Colv = FALSE, col = colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255))
dev.off()


plotPCA(rld_sub)
(data <- plotPCA(rld_sub, returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))
library("ggplot2")
shapes <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,25)
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + scale_shape_manual(values=shapes) + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())

pdf(".../PCA_3D_conditions_shapesize4.pdf")
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + scale_shape_manual(values=shapes) + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
dev.off()

pdf(".../PCA_3D_conditions_shapesize5.pdf")
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=5) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + scale_shape_manual(values=shapes) + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
dev.off()

pdf(".../PCA_3D_conditions_shapesize6.pdf")
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=6) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + scale_shape_manual(values=shapes) + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
dev.off()

pdf(".../PCA_3D_conditions_shapesize7.pdf")
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=7) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + scale_shape_manual(values=shapes) + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
dev.off()
