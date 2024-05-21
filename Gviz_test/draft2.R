install.packages('openxlsx')

library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)

# Assembly
gen <- "hg19"


# MECP2
mecp2_chr  <- "chrX"
mecp2_from <- 153287024
mecp2_to   <- 153363212
mecp2_pos  <- 153363075

# WDR45
wd45_chr  <- "chrX"
wd45_from <- 153287024
wd45_to   <- 153363212
wd45_pos  <- 153286024



# ideogram
itrack <- IdeogramTrack(genome = gen, chromosome = chr, cex = 1, 
                        bevel = 1, showId = FALSE)

# axis
axistrack <- GenomeAxisTrack(add35 = F, add53 = F, cex = 1, 
                             fontcolor = "#383838")

# Genes
ucscTrack1 <- UcscTrack(genome=gen, chromosome=chr, track="NCBI RefSeq", 
                   table = "refGene", from=from, to=to, rstarts="exonStarts", 
                   rends="exonEnds", trackType="GeneRegionTrack", 
                   gene="name", symbol="name", name="Genes",
                   transcript = "name", strand = "strand"
                   )

wdr45 <- UcscTrack(genome=gen, chromosome=chr, track="NCBI RefSeq", 
                   table = "refGene",from=from, to=to, rstarts="exonStarts", 
                   rends="exonEnds", trackType="GeneRegionTrack", 
                   gene="name", symbol="name", name="Genes",
                   transcript = "name", strand = "strand"
                   )

displayPars(mecp2) <- list(fill = "#665990", background.title = "#665990", 
                           stacking = 'full', transcriptAnnotation = "transcript", col = NULL
                           )

tracklist = list(itrack, mecp2)
ht <- HighlightTrack(trackList = tracklist, alpha = 0.4,
                     start = mecp2_pos, width = 2000,
                     chromosome = chr, inBackground = FALSE)

plotTracks(ht, from=from-20000, to=to+20000)

tracklist = list(itrack, mecp2, ht)


zoom in intersetibg




# sequence
strack <- SequenceTrack(Hsapiens, chromosome = chr, 
                        add53 = TRUE, complement = FALSE)
plotTracks(list(ncbirefseq, strack), from = from, to = to)

# axis
gtrack <- GenomeAxisTrack(add35 = TRUE, add53 = TRUE, cex = 1.1,
                          fontcolor = "#383838")

# ideogram
itrack <- IdeogramTrack(genome = gen, chromosome = chr, cex = 0.5, bevel = 1, 
                        showId = FALSE)

plotTracks(list(itrack, gtrack, knownGenes, strack), from = from, to = to)

# Data
data <- data.frame(start = c(159682600, 159682650), end = c(159682600, 159682650), 
                   deltrascore = c(-0.91, -0.65))


dtrack <- DataTrack(range = data, chromosome = chr, genome=gen, 
                    name="SpliceAI delta score", cex.title = 1.0, 
                    background.title = "#F8ACAC", type = "histogram", 
                    baseline = 0, ylim = c(-1, 1), lwd.baseline = 1, 
                    yTicksAt = c(-1.0, -0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8, 1.0))


plotTracks(list(itrack, gtrack, knownGenes, dtrack, strack), 
           from = 159682200, to = 153363000)

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
