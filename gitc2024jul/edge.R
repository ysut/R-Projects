if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
library(edgeR)

dat <- read.delim("/Volumes/vol/gitc/data/SS/arab2.txt", row.names=1)

# assign groups
grp <- c("M", "M", "M", "H", "H", "H")
grp

# import data to edgeR
y <- DGEList(counts=dat, group=grp)
y$samples

# normalization (TMM)
y <- normLibSizes(y, method="TMM")
y$samples

# estimate dispersion
y <- estimateDisp(y)
y$common.dispersion

# DE test
et <- exactTest(y)
et

topTags(et, n=10, adjust.method="BH")
et.sorted <- topTags(et, n=nrow(et$table), adjust.method="BH")
write.table(et.sorted$table, "et.sorted.txt", sep="\t", quote=F)

et.sorted <- topTags(et, n=nrow(et$table))
detab <- et.sorted$table
detab[detab$logFC > log2(10),]

nrow(detab[detab$logFC > log2(10),])



detab[(detab$logFC > log2(5)  & detab$FDR < 0.01), ]

summary(decideTests(et, p.value=0.01, lfc=log2(5)))
plotSmear(y)

de.names <- row.names(et[decideTestsDGE(et, p.value=0.05, lfc=log2(2)) !=0, ])
plotSmear(y, de.tags=de.names)














