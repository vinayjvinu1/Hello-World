
snp <- read.table("SNPs.txt", header = T, row.names = NULL, sep = "\t")
View(snp)

file <- read.table("1T.txt", header = T, sep = "\t")
View(file)


mergedata2 <-merge(file, snp, by.x = "POS", by.y = "POS", all.x = TRUE)
View(mergedata2)

write.csv(mergedata2,"merge2.csv")

m1 <- merge(snp, file, by.x = "POS")
