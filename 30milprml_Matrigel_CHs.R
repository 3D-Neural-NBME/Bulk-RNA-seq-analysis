
countData = read.table(file=".../rsem.genes.counts.matrix",header=TRUE,sep="\t",row.names=1,check.names = FALSE)
colData = read.table(file=".../samples_described.txt",head=TRUE,sep="\t",row.names=1,check.names = FALSE)
rownames(colData)=sub("-","_",rownames(colData), fixed=T)
countData <- countData[ , rownames(colData)]
countData = round(countData)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

dds_sub <- dds[, dds$condition %in% c("ConK", "Con17X", "Con3X", "Con2X", "ConP")]
dds_sub$condition <- droplevels(dds_sub$condition)
rld_sub <- rlog(dds_sub, blind=FALSE)

plotPCA(rld_sub)
(data <- plotPCA(rld_sub, returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))
library("ggplot2")
#shapes <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)
#shapes <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,25)
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())

pdf(".../PCA_ConK_17X_3X_2X_P.pdf")
ggplot(data, aes(PC1, PC2, color=condition, shape=condition))+ geom_point(size=4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed() + theme_bw() +  theme(axis.text = element_text(size = 10)) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
dev.off()