library(biomaRt)
library(openxlsx)
db <- useMart("ensembl")
hd <- useDataset("hsapiens_gene_ensembl", mart = db)
ids = c("hgnc_symbol", "hgnc_id", "description")

setwd("/Volumes/SSD_480GB/workspace/Github/R-Projects/IDconvert")

listDatasets(db)
filters <- listFilters(hd)
filters[grep("prev", filters[,1]),]

# hgnc_symbol
# hgnc_id
# entrezgene_id
# ensembl_gene_id
# ensembl_transcript_id_version
# uniprot_gn_symbol
# refseq_mrna

attributes <- listAttributes(hd)
attributes[grep("prev", attributes[,1]),]

#### TEST SPACE ####
testdata <- c("KMT2B", "KMT2D", "MTTP", "ARRS", "MT-TP")
res <- getBM(attributes = ids,
             filters = "hgnc_symbol", values = testdata, 
             mart = hd, useCache = FALSE)
res
####################

################################################################################
data <- read.xlsx("ajhg.xlsx", 1)
res <- getBM(attributes = ids,
             filters = "hgnc_symbol", values = data$gene, 
             mart = hd, useCache = FALSE)

res
write.table(res, "ajhg.tsv", sep = "\t", row.names = FALSE)

################################################################################
higene_xlsx <- "/Volumes/SSD_480GB/workspace/Github/dev/Resources/02_EstimatedLoFGenes/higene.xlsx"
data <- read.xlsx(higene_xlsx, 1)
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
phaplo_xlsx <- "/Volumes/SSD_480GB/workspace/Github/dev/Resources/02_EstimatedLoFGenes/phaplo.xlsx"
phaplo_xlsx <- "/Volumes/SSD_480GB/workspace/Github/dev/Resources/02_EstimatedLoFGenes/phaplo_chk.xlsx"
data <- read.xlsx(phaplo_xlsx, 1)
res4 <- getBM(attributes = ids,
              filters = "hgnc_symbol", values = data$gene, 
              mart = hd, useCache = FALSE)

res4
write.table(res4, "phaplo.tsv", sep = "\t", row.names = FALSE)

################################################################################
data <- read.xlsx("genovo.xlsx", 1)
ids = c("ensembl_transcript_id", "hgnc_symbol", "hgnc_id", "description")
res5 <- getBM(attributes = ids,
              filters = "ensembl_transcript_id", values = data$enstID, 
              mart = hd, useCache = FALSE)

res5
write.table(res5, "genovo.tsv", sep = "\t", row.names = FALSE)


################################################################################
data <- read.xlsx("gnocchi.xlsx", 1)
ids = c("hgnc_symbol", "hgnc_id", "description")
res6 <- getBM(attributes = ids,
              filters = "hgnc_symbol", values = data$gene, 
              mart = hd, useCache = FALSE)
res6
write.table(res6, "gnocchi.tsv", sep = "\t", row.names = FALSE)

################################################################################
data <- read.xlsx("am.xlsx", 1)
ids = c("hgnc_symbol", "hgnc_id", "description")
res7 <- getBM(attributes = ids,
              filters = "hgnc_symbol", values = data$gene, 
              mart = hd, useCache = FALSE)
res7
write.table(res7, "am.tsv", sep = "\t", row.names = FALSE)


