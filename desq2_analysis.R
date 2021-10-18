library(DESeq2)
library(apeglm)

#load count matrix

dat <- read.table("Gall.txt", header = T, row.names = 1, sep = "\t")
View(dat)

info <- read.table("man.txt", header = T, sep = ",")
View(info)


dds <- DESeqDataSetFromMatrix(dat, info, design=~Tissue, tidy = TRUE)

#remove lowly expressed genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#main DESeq 
ddsDE <- DESeq(dds)

#export norma read counts
normCounts <- counts(ddsDE, normalized = T)
write.csv(normCounts, "normal.csv")

#DESeq results

res <- results(ddsDE, alpha = 0.05)
summary(res)

resOrdered <- res[order(res$padj),]
write.csv(resOrdered,"deSeq.csv")

plotMA(ddsDE, ylim = c(-5,5))

