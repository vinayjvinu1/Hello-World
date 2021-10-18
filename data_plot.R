
library(ggplot2)
library(pheatmap)

normCount <- read.csv("normal.csv", row.names = 1)
View(normCount)

deSeqRes <- read.csv("deSeq.csv", row.names = 1)

deSeqRes$sig <- ifelse(deSeqRes$padj <= 0.05, "yes", "no")

deSeqRes <- na.omit(deSeqRes)

ggplot(deSeqRes, aes(x = log10(baseMean), y = log2FoldChange, color = sig)) +
  geom_point()

#volcano plot

ggplot(deSeqRes, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point()


#pheat map

signi <- subset(deSeqRes, padj <= 0.05)

allsig <- merge(normCount, signi, by = 0)
sigCounts <- allsig[,2:9]

row.names(sigCounts) <- allsig$Row.names

pheatmap(sigCounts)

pheatmap(log2(sigCounts +1))

pheatmap(log2(sigCounts +1), scale = 'row', show_rownames = F)



