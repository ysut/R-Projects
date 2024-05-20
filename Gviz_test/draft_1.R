library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)

install.packages('openxlsx')

chr <- "chr1"
gen <- "hg19"
from <- 159681503
to   <- 159684954

knownGenes <- UcscTrack(genome=gen, chromosome=chr, track="NCBI RefSeq", 
                        table = "ncbiRefSeqSelect",
                        from=from, to=to,trackType="GeneRegionTrack", 
                        rstarts="exonStarts", rends="exonEnds", 
                        gene="name", symbol="name", 
                        transcript="name", strand="strand", 
                        fill="#AAA5D1", name="UCSC Genes")
displayPars(knownGenes) <- list(background.title = "#665990")

# sequence
strack <- SequenceTrack(Hsapiens, chromosome = chr, 
                        add53 = TRUE, complement = FALSE)

# axis
gtrack <- GenomeAxisTrack(add35 = TRUE, add53 = TRUE, cex = 1.1,
                          fontcolor = "#383838")

# ideogram
itrack <- IdeogramTrack(genome = gen, chromosome = chr, cex = 0.5, bevel = 1, 
                        showId = FALSE)

# Data
data <- data.frame(start = c(159682600, 159682650), end = c(159682600, 159682650), 
                   deltrascore = c(-0.91, -0.65))


dtrack <- DataTrack(range = data, chromosome = chr, genome=gen, 
                    name="SpliceAI delta score", cex.title = 1.0, 
                    background.title = "#F8ACAC", type = "histogram", 
                    baseline = 0, ylim = c(-1, 1), lwd.baseline = 1, 
                    yTicksAt = c(-1.0, -0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8, 1.0))

plotTracks(list(itrack, gtrack, knownGenes, dtrack, strack), from = 159682200, to = 159682680)

library(Gviz)
library(openxlsx)
gene_model_xlsx <- '/Users/utsu/Google Drive/マイドライブ/01_Manuscripts/SplicingScreening/99_Figure_drafts/gene_model.xlsx'
scores_xlsx <- '/Users/utsu/Google Drive/マイドライブ/01_Manuscripts/SplicingScreening/99_Figure_drafts/scores.xlsx'

atrack <- AnnotationTrack(start = 270, width = 0, name = 'Variant Pos.',
                          group = "Variant", background.title = "#383838", shape = "box", 
                          genome = "hg19", chromosome = "chr1", just.group = "right", showOverplotting = TRUE)

axtrack <- GenomeAxisTrack(add35 = TRUE, add53 = TRUE, cex = 1,
                           fontcolor = "#383838")

scores <- read.xlsx(scores_xlsx)
splaitrack <- DataTrack(range = scores, chromosome = "chr1", genome="hg19", 
                        name="SpliceAI delta score", cex.title = 1.0, 
                        background.title = "#F8ACAC", type = "histogram", 
                        baseline = 0, ylim = c(-1, 1), lwd.baseline = 1, 
                        yTicksAt = c(-1.0, -0.5, -0.2, 0, 0.2, 0.5, 1.0),
                        groups = c("AG", "AL", "DG", "DL"),
                        col = c("#6088C6", "#EB8686", "#73D0C2", "#ED8D49"),
                        cex.legend = 1.0, box.legend = FALSE)


genemodel <- read.xlsx(gene_model_xlsx)
gmtrack <- GeneRegionTrack(genemodel, name = "Gene Model", background.title = "#665990", fill = "#665990",
                           chromosome = "chr1", genome = "hg19")

# displayPars(gmtrack) <- list(alpha.title = 1, alpha = 0.5)
plotTracks(list(axtrack, atrack, gmtrack, splaitrack), from = 0, to = 600)

ot <- OverlayTrack(trackList=list(gmtrack, atrack))
plotTracks(ot, from = 0, to = 700)
