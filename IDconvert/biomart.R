library(biomaRt)
library(openxlsx)
db <- useMart("ensembl")
hd <- useDataset("hsapiens_gene_ensembl", mart = db)
ids = c("hgnc_symbol", "hgnc_id", "entrezgene_id", "description")


listDatasets(db)
filters <- listFilters(hd)
filters[grep("transcript", filters[,1]),]

# hgnc_symbol
# hgnc_id
# entrezgene_id
# ensembl_gene_id
# ensembl_transcript_id_version
# uniprot_gn_symbol
# refseq_mrna


################################################################################
data <- read.xlsx("ajhg.ori.xlsx", 1)
res <- getBM(attributes = ids,
             filters = "hgnc_symbol", values = data$gene, 
             mart = hd, useCache = FALSE)

res
write.table(res, "ajhg.tsv", sep = "\t", row.names = FALSE)

################################################################################
data <- read.xlsx("higene.xlsx", 1)
res2 <- getBM(attributes = ids,
              filters = "hgnc_symbol", values = data$gene, 
              mart = hd, useCache = FALSE)
res2
write.table(res2, "higene.tsv", sep = "\t", row.names = FALSE)

################################################################################
data <- read.xlsx("gnomad.xlsx", 1)
res3 <- getBM(attributes = ids,
              filters = "hgnc_symbol", values = data$gene, 
              mart = hd, useCache = FALSE)
res3
write.table(res3, "gnomad.tsv", sep = "\t", row.names = FALSE)

################################################################################
data <- read.xlsx("phaplo.xlsx", 1)
res4 <- getBM(attributes = ids,
              filters = "hgnc_symbol", values = data$gene, 
              mart = hd, useCache = FALSE)

res4
write.table(res4, "phaplo.tsv", sep = "\t", row.names = FALSE)

################################################################################
data <- read.xlsx("genovo.xlsx", 1)
ids = c("ensembl_transcript_id", "hgnc_symbol", "hgnc_id", "entrezgene_id", "description")
res5 <- getBM(attributes = ids,
              filters = "ensembl_transcript_id", values = data$enstID, 
              mart = hd, useCache = FALSE)

res5
write.table(res5, "genovo.tsv", sep = "\t", row.names = FALSE)


################################################################################
data <- read.xlsx("gnocchi.xlsx", 1)
ids = c("hgnc_symbol", "hgnc_id", "entrezgene_id", "description")
res6 <- getBM(attributes = ids,
              filters = "hgnc_symbol", values = data$gene, 
              mart = hd, useCache = FALSE)
res6
write.table(res6, "gnocchi.tsv", sep = "\t", row.names = FALSE)

################################################################################
data <- read.xlsx("am.xlsx", 1)
ids = c("hgnc_symbol", "hgnc_id", "entrezgene_id", "description")
res7 <- getBM(attributes = ids,
              filters = "hgnc_symbol", values = data$gene, 
              mart = hd, useCache = FALSE)
res7
write.table(res7, "am.tsv", sep = "\t", row.names = FALSE)


