library(biomaRt)
library(openxlsx)
db <- useMart("ensembl")
hd <- useDataset("hsapiens_gene_ensembl", mart = db)
ids = c("hgnc_symbol", "hgnc_id", "entrezgene_id", "description")


listDatasets(db)
filters <- listFilters(hd)
filters[grep("symbol", filters[,1]),]

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

################################################################################












