
countData = read.table(file=".../rsem.genes.counts.matrix",header=TRUE,sep="\t",row.names=1,check.names = FALSE)
colData = read.table(file=".../samples_described.txt",head=TRUE,sep="\t",row.names=1,check.names = FALSE)
rownames(colData)=sub("-","_",rownames(colData), fixed=T)
countData <- countData[ , rownames(colData)]
countData = round(countData)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)


dds_sub <- dds[, dds$condition %in% c("Con1Z","Con5Z","Con6Z","Con7Z","Con8Z","Con9Z","Con10Z","Con2Z","Con13Z","Con11Z","Con12Z","Con14Z","Con15Z","Con3Z")]
dds_sub$condition <- droplevels(dds_sub$condition)
rld_sub <- rlog(dds_sub, blind=FALSE)


## for supplementary figure
## z-score heatmap with Inhibitory neuronal markers, Neural epithelial markers genes

selected_genes <- c("GFAP","S100B","ALDH1L1","AQP4","AGT","SOX9","VIM","SLC1A3","PAX6","SOX2","HES1","GLI3","KLF4", "VCAN", "COL1A2","EOMES","PPP1R17","PENK","NEUROD1","HES6","ELAVL4","SATB2","SLC17A7","SLC17A6","GRIN1","GRIN2A","GRIN2B", "GAD1", "GAD2", "SLC32A1", "DLX2", "DLX5", "DLX6")

mat <- assay(rld_sub)[selected_genes,]

mat <- mat[,c("1A_070517_S1","1B_070517_S2","1C_070517_S3","5A_070517_S13","5B_070517_S14","5C_070517_S15","6A_070517_S16","6B_070517_S17","6C_070517_S18","7A_070517_S19","7B_070517_S20","7C_070517_S21","8A_070517_S22","8B_070517_S23","8C_070517_S24",
              "9A_070517_S25","9B_070517_S26","9C_070517_S27","10A_070517_S28","10B_070517_S29","10C_070517_S30","2A_070517_S4","2B_070517_S5","2C_070517_S6","13A_070517_S37","13B_070517_S38","13C_070517_S39","11A_070517_S31","11B_070517_S32","11C_070517_S34",
              "12A_070517_S33","12B_070517_S35","12C_070517_S36","14A_070517_S40","14B_070517_S41","14C_070517_S42","15A_070517_S43","15B_070517_S44","15C_070517_S45","3A_070517_S7","3B_070517_S8","3C_070517_S9")]

library("gplots")
library("RColorBrewer")
heatmap.2(mat,scale="row",trace="none", dendrogram="none", Rowv = FALSE,Colv = FALSE, col = colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255))

pdf(".../for_suppfigure_selected_marker_genes_with_zscore.pdf")
heatmap.2(mat,scale="row",trace="none", dendrogram="none", Rowv = FALSE,Colv = FALSE, col = colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255))
dev.off()



### for main figure - subset of conditions hESC, day2, day30, day67, day114, huPAst ##
## z-score heatmap with Inhibitory neuronal markers, Neural epithelial markers genes

dds_sub <- dds[, dds$condition %in% c("Con1Z","Con5Z","Con13Z","Con11Z","Con12Z","Con14Z","Con15Z","Con3Z")]
dds_sub$condition <- droplevels(dds_sub$condition)
rld_sub <- rlog(dds_sub, blind=FALSE)

selected_genes <- c("GFAP","S100B","ALDH1L1","AQP4","AGT","SOX9","VIM","SLC1A3","PAX6","SOX2","HES1","GLI3","KLF4", "VCAN", "COL1A2","EOMES","PPP1R17","PENK","NEUROD1","HES6","ELAVL4","SATB2","SLC17A7","SLC17A6","GRIN1","GRIN2A","GRIN2B", "GAD1", "GAD2", "SLC32A1", "DLX2", "DLX5", "DLX6")

mat <- assay(rld_sub)[selected_genes,]

mat <- mat[,c("1A_070517_S1","1B_070517_S2","1C_070517_S3","5A_070517_S13","5B_070517_S14","5C_070517_S15", "13A_070517_S37","13B_070517_S38","13C_070517_S39","11A_070517_S31","11B_070517_S32","11C_070517_S34",
              "12A_070517_S33","12B_070517_S35","12C_070517_S36","14A_070517_S40","14B_070517_S41","14C_070517_S42","15A_070517_S43","15B_070517_S44","15C_070517_S45","3A_070517_S7","3B_070517_S8","3C_070517_S9")]

library("gplots")
library("RColorBrewer")
heatmap.2(mat,scale="row",trace="none", dendrogram="none", Rowv = FALSE,Colv = FALSE, col = colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255))

pdf(".../for_mainfigure_subsetconditions_selected_marker_genes_zscore.pdf")
heatmap.2(mat,scale="row",trace="none", dendrogram="none", Rowv = FALSE,Colv = FALSE, col = colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255))
dev.off()


