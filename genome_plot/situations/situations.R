library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(openxlsx)

# Genome assembly
gen <- "hg19"

gene_model_xlsx <- "situations.xlsx"


genemodel <- read.xlsx(gene_model_xlsx)
gmtrack <- GeneRegionTrack(genemodel, name = "Gene Model", background.title = "#665990", fill = "#665990",
                           chromosome = "chr1", genome = "hg19")