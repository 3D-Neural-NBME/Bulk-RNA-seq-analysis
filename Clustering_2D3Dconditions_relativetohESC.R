
countData = read.table(file=".../rsem.genes.counts.matrix",header=TRUE,sep="\t",row.names=1,check.names = FALSE)
colData = read.table(file=".../samples_described.txt",head=TRUE,sep="\t",row.names=1,check.names = FALSE)
rownames(colData)=sub("-","_",rownames(colData), fixed=T)
countData <- countData[ , rownames(colData)]
countData = round(countData)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

dds_sub <- dds[, dds$condition %in% c("ConL","Con20X","ConN","ConI","ConG","Con18Z","ConJ","ConH","ConK","ConO","Con18X","Con17X","Con3X","Con2X","ConP")]

dds_Con20X <- dds[, dds$condition %in% c("ConL", "Con20X")]
dds_Con20X$condition <- relevel(dds_Con20X$condition, ref="ConL")
dds_Con20X$condition <- droplevels(dds_Con20X$condition)
dds_Con20X <- DESeq(dds_Con20X)
res05_Con20X <- results(dds_Con20X, alpha=0.05)
res05Ordered_Con20X <- res05_Con20X[order(res05_Con20X$padj),]
data_Con20X <- na.omit(res05Ordered_Con20X)
fltr_Con20X <- data_Con20X[(as.matrix(data_Con20X[,6]) < 0.01)&(as.matrix(data_Con20X[,2]) > 2 | as.matrix(data_Con20X[,2]) < -2),]
dim(fltr_Con20X)

dds_ConN <- dds[, dds$condition %in% c("ConL", "ConN")]
dds_ConN$condition <- relevel(dds_ConN$condition, ref="ConL")
dds_ConN$condition <- droplevels(dds_ConN$condition)
dds_ConN <- DESeq(dds_ConN)
res05_ConN <- results(dds_ConN, alpha=0.05)
res05Ordered_ConN <- res05_ConN[order(res05_ConN$padj),]
data_ConN <- na.omit(res05Ordered_ConN)
fltr_ConN <- data_ConN[(as.matrix(data_ConN[,6]) < 0.01)&(as.matrix(data_ConN[,2]) > 2 | as.matrix(data_ConN[,2]) < -2),]
dim(fltr_ConN)

dds_ConI <- dds[, dds$condition %in% c("ConL", "ConI")]
dds_ConI$condition <- relevel(dds_ConI$condition, ref="ConL")
dds_ConI$condition <- droplevels(dds_ConI$condition)
dds_ConI <- DESeq(dds_ConI)
res05_ConI <- results(dds_ConI, alpha=0.05)
res05Ordered_ConI <- res05_ConI[order(res05_ConI$padj),]
data_ConI <- na.omit(res05Ordered_ConI)
fltr_ConI <- data_ConI[(as.matrix(data_ConI[,6]) < 0.01)&(as.matrix(data_ConI[,2]) > 2 | as.matrix(data_ConI[,2]) < -2),]
dim(fltr_ConI)

dds_ConG <- dds[, dds$condition %in% c("ConL", "ConG")]
dds_ConG$condition <- relevel(dds_ConG$condition, ref="ConL")
dds_ConG$condition <- droplevels(dds_ConG$condition)
dds_ConG <- DESeq(dds_ConG)
res05_ConG <- results(dds_ConG, alpha=0.05)
res05Ordered_ConG <- res05_ConG[order(res05_ConG$padj),]
data_ConG <- na.omit(res05Ordered_ConG)
fltr_ConG <- data_ConG[(as.matrix(data_ConG[,6]) < 0.01)&(as.matrix(data_ConG[,2]) > 2 | as.matrix(data_ConG[,2]) < -2),]
dim(fltr_ConG)

dds_Con18Z <- dds[, dds$condition %in% c("ConL", "Con18Z")]
dds_Con18Z$condition <- relevel(dds_Con18Z$condition, ref="ConL")
dds_Con18Z$condition <- droplevels(dds_Con18Z$condition)
dds_Con18Z <- DESeq(dds_Con18Z)
res05_Con18Z <- results(dds_Con18Z, alpha=0.05)
res05Ordered_Con18Z <- res05_Con18Z[order(res05_Con18Z$padj),]
data_Con18Z <- na.omit(res05Ordered_Con18Z)
fltr_Con18Z <- data_Con18Z[(as.matrix(data_Con18Z[,6]) < 0.01)&(as.matrix(data_Con18Z[,2]) > 2 | as.matrix(data_Con18Z[,2]) < -2),]
dim(fltr_Con18Z)

dds_ConJ <- dds[, dds$condition %in% c("ConL", "ConJ")]
dds_ConJ$condition <- relevel(dds_ConJ$condition, ref="ConL")
dds_ConJ$condition <- droplevels(dds_ConJ$condition)
dds_ConJ <- DESeq(dds_ConJ)
res05_ConJ <- results(dds_ConJ, alpha=0.05)
res05Ordered_ConJ <- res05_ConJ[order(res05_ConJ$padj),]
data_ConJ <- na.omit(res05Ordered_ConJ)
fltr_ConJ <- data_ConJ[(as.matrix(data_ConJ[,6]) < 0.01)&(as.matrix(data_ConJ[,2]) > 2 | as.matrix(data_ConJ[,2]) < -2),]
dim(fltr_ConJ)

dds_ConH <- dds[, dds$condition %in% c("ConL", "ConH")]
dds_ConH$condition <- relevel(dds_ConH$condition, ref="ConL")
dds_ConH$condition <- droplevels(dds_ConH$condition)
dds_ConH <- DESeq(dds_ConH)
res05_ConH <- results(dds_ConH, alpha=0.05)
res05Ordered_ConH <- res05_ConH[order(res05_ConH$padj),]
data_ConH <- na.omit(res05Ordered_ConH)
fltr_ConH <- data_ConH[(as.matrix(data_ConH[,6]) < 0.01)&(as.matrix(data_ConH[,2]) > 2 | as.matrix(data_ConH[,2]) < -2),]
dim(fltr_ConH)

dds_ConK <- dds[, dds$condition %in% c("ConL", "ConK")]
dds_ConK$condition <- relevel(dds_ConK$condition, ref="ConL")
dds_ConK$condition <- droplevels(dds_ConK$condition)
dds_ConK <- DESeq(dds_ConK)
res05_ConK <- results(dds_ConK, alpha=0.05)
res05Ordered_ConK <- res05_ConK[order(res05_ConK$padj),]
data_ConK <- na.omit(res05Ordered_ConK)
fltr_ConK <- data_ConK[(as.matrix(data_ConK[,6]) < 0.01)&(as.matrix(data_ConK[,2]) > 2 | as.matrix(data_ConK[,2]) < -2),]
dim(fltr_ConK)

dds_ConO <- dds[, dds$condition %in% c("ConL", "ConO")]
dds_ConO$condition <- relevel(dds_ConO$condition, ref="ConL")
dds_ConO$condition <- droplevels(dds_ConO$condition)
dds_ConO <- DESeq(dds_ConO)
res05_ConO <- results(dds_ConO, alpha=0.05)
res05Ordered_ConO <- res05_ConO[order(res05_ConO$padj),]
data_ConO <- na.omit(res05Ordered_ConO)
fltr_ConO <- data_ConO[(as.matrix(data_ConO[,6]) < 0.01)&(as.matrix(data_ConO[,2]) > 2 | as.matrix(data_ConO[,2]) < -2),]
dim(fltr_ConO)

dds_Con18X <- dds[, dds$condition %in% c("ConL", "Con18X")]
dds_Con18X$condition <- relevel(dds_Con18X$condition, ref="ConL")
dds_Con18X$condition <- droplevels(dds_Con18X$condition)
dds_Con18X <- DESeq(dds_Con18X)
res05_Con18X <- results(dds_Con18X, alpha=0.05)
res05Ordered_Con18X <- res05_Con18X[order(res05_Con18X$padj),]
data_Con18X <- na.omit(res05Ordered_Con18X)
fltr_Con18X <- data_Con18X[(as.matrix(data_Con18X[,6]) < 0.01)&(as.matrix(data_Con18X[,2]) > 2 | as.matrix(data_Con18X[,2]) < -2),]
dim(fltr_Con18X)

dds_Con17X <- dds[, dds$condition %in% c("ConL", "Con17X")]
dds_Con17X$condition <- relevel(dds_Con17X$condition, ref="ConL")
dds_Con17X$condition <- droplevels(dds_Con17X$condition)
dds_Con17X <- DESeq(dds_Con17X)
res05_Con17X <- results(dds_Con17X, alpha=0.05)
res05Ordered_Con17X <- res05_Con17X[order(res05_Con17X$padj),]
data_Con17X <- na.omit(res05Ordered_Con17X)
fltr_Con17X <- data_Con17X[(as.matrix(data_Con17X[,6]) < 0.01)&(as.matrix(data_Con17X[,2]) > 2 | as.matrix(data_Con17X[,2]) < -2),]
dim(fltr_Con17X)

dds_Con3X <- dds[, dds$condition %in% c("ConL", "Con3X")]
dds_Con3X$condition <- relevel(dds_Con3X$condition, ref="ConL")
dds_Con3X$condition <- droplevels(dds_Con3X$condition)
dds_Con3X <- DESeq(dds_Con3X)
res05_Con3X <- results(dds_Con3X, alpha=0.05)
res05Ordered_Con3X <- res05_Con3X[order(res05_Con3X$padj),]
data_Con3X <- na.omit(res05Ordered_Con3X)
fltr_Con3X <- data_Con3X[(as.matrix(data_Con3X[,6]) < 0.01)&(as.matrix(data_Con3X[,2]) > 2 | as.matrix(data_Con3X[,2]) < -2),]
dim(fltr_Con3X)

dds_Con2X <- dds[, dds$condition %in% c("ConL", "Con2X")]
dds_Con2X$condition <- relevel(dds_Con2X$condition, ref="ConL")
dds_Con2X$condition <- droplevels(dds_Con2X$condition)
dds_Con2X <- DESeq(dds_Con2X)
res05_Con2X <- results(dds_Con2X, alpha=0.05)
res05Ordered_Con2X <- res05_Con2X[order(res05_Con2X$padj),]
data_Con2X <- na.omit(res05Ordered_Con2X)
fltr_Con2X <- data_Con2X[(as.matrix(data_Con2X[,6]) < 0.01)&(as.matrix(data_Con2X[,2]) > 2 | as.matrix(data_Con2X[,2]) < -2),]
dim(fltr_Con2X)

dds_ConP <- dds[, dds$condition %in% c("ConL", "ConP")]
dds_ConP$condition <- relevel(dds_ConP$condition, ref="ConL")
dds_ConP$condition <- droplevels(dds_ConP$condition)
dds_ConP <- DESeq(dds_ConP)
res05_ConP <- results(dds_ConP, alpha=0.05)
res05Ordered_ConP <- res05_ConP[order(res05_ConP$padj),]
data_ConP <- na.omit(res05Ordered_ConP)
fltr_ConP <- data_ConP[(as.matrix(data_ConP[,6]) < 0.01)&(as.matrix(data_ConP[,2]) > 2 | as.matrix(data_ConP[,2]) < -2),]
dim(fltr_ConP)

rows_20X <- rownames(fltr_Con20X)
rows_N <- rownames(fltr_ConN)
rows_I <- rownames(fltr_ConI)
rows_G <- rownames(fltr_ConG)
rows_18Z <- rownames(fltr_Con18Z)
rows_J <- rownames(fltr_ConJ)
rows_H <- rownames(fltr_ConH)
rows_K <- rownames(fltr_ConK)
rows_O <- rownames(fltr_ConO)
rows_18X <- rownames(fltr_Con18X)
rows_17X <- rownames(fltr_Con17X)
rows_3X <- rownames(fltr_Con3X)
rows_2X <- rownames(fltr_Con2X)
rows_P <- rownames(fltr_ConP)

comb_20XN <- c(rows_20X,rows_N)
comb_IG <- c(rows_I,rows_G)
comb_18ZJ <- c(rows_18Z,rows_J)
comb_HK <- c(rows_H,rows_K)
comb_O18X <- c(rows_O,rows_18X)
comb_17X3X <- c(rows_17X,rows_3X)
comb_2XP <- c(rows_2X,rows_P)

comb_a <- c(comb_20XN,comb_IG)
comb_b <- c(comb_18ZJ,comb_HK)
comb_c <- c(comb_O18X,comb_17X3X)
comb_ab <- c(comb_a,comb_b)
comb_abc <- c(comb_ab,comb_c)
combine_all <- c(comb_2XP,comb_abc)
genes_all <- unique(combine_all)
length(genes_all)

write.csv(as.data.frame(genes_all), file=".../DErelativetohESC_14conditions_combineduniquegenelist_p0.01_log2fc2_orminus2.csv")

dds_sub$condition <- droplevels(dds_sub$condition)
rld_sub <- rlog(dds_sub, blind=FALSE)
mat <- assay(rld_sub)[genes_all,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds_sub)[c("condition")])
library("pheatmap")
library("RColorBrewer")

km<- kmeans(mat,4)
m.kmeans<- cbind(mat, km$cluster)
dim(m.kmeans)
o<- order(m.kmeans[,45])
m.kmeans<- m.kmeans[o,]

a<-m.kmeans[,1:44]
a <- a[,c("12A_S34","12B_S35","12C_S36","20A061516_S61","20B061516_S62","14A_S40","14B_S41","14C_S42","9A_S25","9B_S26","9C_S27","7A_S19","7B_S20","7C_S21","18A_070517_S52","18B_070517_S53","18C_070517_S54",
          "10A_S28","10B_S29","10C_S30", "8A_S22","8B_S23","8C_S24","11A_S31","11B_S32","11C_S33","15A_S43","15B_S44","15C_S45","18A060216_S51","18B060216_S52","18C060216_S53",
          "17A060216_S48","17B060216_S49","17C060216_S50","3A060116_S7","3B060116_S8","3C060116_S9","2A060116_S4","2B060116_S5","2C060116_S6","16A_S46","16B_S47","16C_S48")]
          
pheatmap(a, cluster_rows = F, cluster_cols = F,  annotation_col=df, show_rownames=FALSE, color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
pdf(".../kmeans_4clusters.pdf")
pheatmap(a, cluster_rows = F, cluster_cols = F,  annotation_col=df, show_rownames=FALSE, color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100))
pheatmap(a, cluster_rows = F, cluster_cols = F,  annotation_col=df, show_rownames=FALSE, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))
dev.off()
write.csv(as.data.frame(a), file=".../kmeans_4clusters_columnsorted.csv")
write.csv(as.data.frame(m.kmeans), file=".../kmeans_4clusters_withclusternumbers.csv")


