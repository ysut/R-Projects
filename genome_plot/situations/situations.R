library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(openxlsx)

# Genome assembly
gen <- "hg19"

gene_model_xlsx <- "situations.xlsx"


genemodel <- read.xlsx(gene_model_xlsx)
gmtrack <- GeneRegionTrack(genemodel, name = "Gene Model", background.title = "#665990", fill = "#665990",
                           chromosome = "chr1", genome = "hg19", size = 1)

axTrack <- GenomeAxisTrack(
  add35 = FALSE, add53 = FALSE, exponent = 0, fontcolor = "#383838",
  fontsize = 10, labelPos = "above", size = 3, distFromAxis = 1.2
)
plotTracks(c(axTrack, gmtrack), from = 960, to = 1300)
plotTracks(c(gmtrack), from = 960, to = 1300)

sTrack <- SequenceTrack(Hsapiens, noLetters=FALSE, chromosome = "chr1") 

plotTracks(c(gmtrack, sTrack), from = 13900, to = 14000, chromosome = "chr1")
