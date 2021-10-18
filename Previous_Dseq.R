library(DESeq2)
library(apeglm)
library(ggplot2)

dat <- read.csv("Gall.txt", header = T, sep = "\t", row.names = 1 )
info <- read.table("man.txt", header = T, sep = ",")

dds <- DESeqDataSetFromMatrix(dat, info, -Tissue)

countData <- read.csv("gall.txt", header = T, sep = "\t")
head(countData)
metaData <- read.csv("man.txt", header = T, sep = ",")
metadata
#Construct DESEQDataSet Object

nams <- read.csv("gall.txt", col.names = 1)


rownames(dat) = make.names(nams, unique=TRUE)

dds <- DESeqDataSetFromMatrix(countData=countData,colData=metaData, 
        design=~Tissue, tidy = TRUE)

dds

#RUN DESEQ FUNCTION
dds <- DESeq(dds)


#Check result table:
  res <- results(dds)
head(results(dds, tidy=TRUE))

summary(res)

#Sort summary list with p-value
res <-res[order(res$padj),]
head(res)

plotCounts
par(mfrow=c(2,3))

plotCounts(dds, gene="ENSG00000152583", intgroup="Tissue")
plotCounts(dds, gene="ENSG00000179094", intgroup="Tissue")
plotCounts(dds, gene="ENSG00000116584", intgroup="Tissue")
plotCounts(dds, gene="ENSG00000189221", intgroup="Tissue")
plotCounts(dds, gene="ENSG00000120129", intgroup="Tissue")
plotCounts(dds, gene="ENSG00000148175", intgroup="Tissue")


#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


all(rownames(coldata) %in% colnames(cts))

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Tissue")])

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
rownames(annotation_c) <- colnames(DAT)
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Tissue")])

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("vsn")



resultsNames(dds)






