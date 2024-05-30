library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg38)
library(openxlsx)
library(txdbmaker)

# Conncecting options
options(Gviz.ucscUrl="http://genome-asia.ucsc.edu/cgi-bin/")

# Genome assembly
gen <- "hg38"

# Case 1
mecp2_chr  <- "chrX"
mecp2_from <- 154021573
mecp2_to   <- 154097717
mecp2_pos <- 154097618
mecp2_bam <- "bams/Sample_11467.X_154021573-154097717.bam"
mecp2_rnaseq <- "bams/142784Aligned.out.sort.MECP2.chrX_154021573_154097717.bam"
mecp2_xlsx <- "splai_mecp2.xlsx"

mecp2_from19 <- 153287024
mecp2_to19   <- 153363174
mecp2_pos19 <- 153363075

# Case 1
chr  <- mecp2_chr
from <- mecp2_from
to   <- mecp2_to
pos  <- mecp2_pos
input_bam <- mecp2_bam
input_rnaseq <- mecp2_rnaseq
scores_xlsx <- mecp2_xlsx

##################################################
# For wide view
##################################################
# ideogram
iTrack <- IdeogramTrack(
  genome = gen, chromosome = chr, cex = 1.2, bevel = 1, showId = TRUE
)

# axis
axTrack <- GenomeAxisTrack(
  add35 = FALSE, add53 = FALSE, exponent = 0, fontcolor = "#383838",
  fontsize = 16, labelPos = "below", size = 3
)

ucscTrack <- UcscTrack(
  genome = gen, track = "NCBI RefSeq", trackType = "GeneRegionTrack",
  rstarts = "exonStarts", rends = "exonEnds", gene ="name", symbol = "name2", 
  transcript = "name", strand = "strand", table = "ncbiRefSeq",
  chromosome = chr, from = from, to = to, name = "Gene"
)

displayPars(ucscTrack) <- list(
  background.title = "#665990", fill ="#665990", cex = 1, 
  transcriptAnnotation = "transcript", fontsize = 12, size = 12
)

for (i in 1:length(start(ucscTrack))) {
  start(ucscTrack)[i] <- start(ucscTrack)[i] + 1
}

# Highlight the variant position
ht <- HighlightTrack(ucscTrack, alpha = 0.5, inBackground = FALSE, 
  chromosome = chr, start = pos - 100, width = 200
)

# Wide view with highlighted variant position
plotTracks(
  list(iTrack, ht, axTrack), from = from - 10000, to = to + 2000, type = "none"
)

plotTracks(
  list(iTrack, ucscTrack, sashimiTrack, axTrack), from = from - 10000, to = to + 2000, type = "none"
)

##################################################
# For wide view with sashimi
##################################################

# Axis 
axTrack <- GenomeAxisTrack(
  add35 = FALSE, add53 = FALSE, exponent = 0, fontcolor = "#383838",
  fontsize = 16, labelPos = "above", size = 8, distFromAxis = 1.2
)
plotTracks(list(axTrack), from = pos - 10, to = pos + 10)


iTrack <- IdeogramTrack(
  genome = gen, chromosome = chr, cex = 2, bevel = 1, showId = TRUE, size = 20
)

#
ucscTrack <- UcscTrack(
  genome = gen, track = "NCBI RefSeq", trackType = "GeneRegionTrack",
  rstarts = "exonStarts", rends = "exonEnds", gene ="name", symbol = "name2", 
  transcript = "name", strand = "strand", table = "ncbiRefSeqCurated",
  chromosome = chr, from = from, to = to, name = "Gene"
)
for (i in 1:length(start(ucscTrack))) {
  start(ucscTrack)[i] <- start(ucscTrack)[i] + 1
}

transcript(ucscTrack) <- c(
  "NM_001369394.2", "NM_001369394.2", "NM_001369394.2", "NM_001316337.2", "NM_001316337.2", "NM_001316337.2", "NM_001316337.2", "NM_001316337.2", "NM_001369391.2", 
  "NM_001369391.2", "NM_001369391.2", "NM_001369391.2", "NM_001369391.2", "NM_001369391.2", "NM_001369391.2", "NM_001369392.2", "NM_001369392.2", "NM_001369392.2", 
  "NM_001369392.2", "NM_001369392.2", "NM_001386137.1", "NM_001386137.1", "NM_001386137.1", "NM_001386137.1", "NM_001386137.1", "NM_001386137.1", "NM_004992.4 (MANE Plus Clinical)", 
  "NM_004992.4 (MANE Plus Clinical)", "NM_004992.4 (MANE Plus Clinical)", "NM_004992.4 (MANE Plus Clinical)", "NM_001386138.1", "NM_001386138.1", "NM_001386138.1", "NM_001386138.1", "NM_001386138.1", "NM_001369393.2", 
  "NM_001369393.2", "NM_001369393.2", "NM_001369393.2", "NM_001110792.2 (MANE Select)", "NM_001110792.2 (MANE Select)", "NM_001110792.2 (MANE Select)", "NM_001386139.1", "NM_001386139.1", "NM_001386139.1", 
  "NM_001386139.1"
)

displayPars(ucscTrack) <- list(
  background.title = "#665990", fill ="#665990", cex = 1.0,
  transcriptAnnotation = "transcript", size = 200, cex.group = 1.3,
  fontsize = 16
)

sashimiTrack <- AlignmentsTrack(
  input_rnaseq, isPaired = TRUE, genome = gen, chromosome = chr, 
  background.title = "darkgrey", type = c("coverage", "sashimi"), 
  stacking = "squish", minCoverageHeight = 2, cex.title = 1.2, fontsize = 16,
  name = "RNA-seq data", coverageHeight = 0.01, size = 1,
  col.axis = "darkgray", sashimiScore = 8, cex.axis = 0
)

# Highlight the variant position
ht <- HighlightTrack(c(ucscTrack, sashimiTrack),alpha = 0.5, inBackground = FALSE, 
                     chromosome = chr, start = pos - 5800, width = 5900
)
plotTracks(c(iTrack, ht), to = to + 1000, from = from - 10000)

##########################
##### Zoomed in view #####
##########################
zoomucscTrack <- UcscTrack(
  genome = gen, track = "NCBI RefSeq", trackType = "GeneRegionTrack",
  rstarts = "exonStarts", rends = "exonEnds", gene ="name", symbol = "name2", 
  transcript = "name", strand = "strand", table = "ncbiRefSeqCurated",
  chromosome = chr, from = pos-500, to = pos+100, name = "Gene"
)

for (i in 1:length(start(zoomucscTrack))) {
  start(zoomucscTrack)[i] <- start(zoomucscTrack)[i] + 1
}

displayPars(zoomucscTrack) <- list(
  background.title = "#665990", fill ="#665990", cex = 1.0, showTitle = FALSE,
  size = 8
)

##
fontcolor = c(A = "#49A190", C = "#0072B2", G = "#1E181B", T = "red", N = 1)
sTrack <- SequenceTrack(
  Hsapiens, size = 10, fontsize = 16, chromosome = chr, fontcolor = fontcolor) 

##
scores <- read.xlsx(scores_xlsx)
splTrack <- DataTrack(
  range = scores, chromosome = chr, genome = gen, name="∆ score", 
  cex.title = 1.0, background.title = "#F8ACAC", type = "histogram", 
  baseline = 0, col.baseline = "#838383", lty.baseline = 2, lwd.baseline = 2,
  ylim = c(-1, 1), lwd.baseline = 2, fontsize = 16,
  yTicksAt = c(-1.0, -0.5,0 , 0.5, 1.0), groups = c("AG", "AL", "DG", "DL"),
  col = c("#6088C6", "#EF8875", "#49A190", "#ED8D49"), legend = TRUE, 
  cex.legend = 1.6, size = 100, cex.axis = 0.9
)

##
alTrack <- AlignmentsTrack(
  input_bam, isPaired = TRUE, genome = gen, chromosome = chr, 
  background.title = "darkgrey", type = c("pileup"), cex.title = 1.0,
  cex = 1, col.axis = "white",
  stacking = "squish", name = "Read data", size = 12, 
  fontsize = 16, cex.axis = 0, 
  max.height = 10, min.height = 0.5, lwd.mismatch = 0.1, alpha.reads = 0.75
)

##
zoomcoverTrack <- AlignmentsTrack(
  input_rnaseq, isPaired = TRUE, genome = gen, chromosome = chr, 
  background.title = "darkgrey", type = c("coverage"), 
  stacking = "squish", minCoverageHeight = 50, cex.title = 1.1, fontsize = 12,
  name = "RNA-seq coverage", coverageHeight = 0.01, size = 5,
  col.axis = "white",
  cex.axis = 0.5
)

## Plot ##
plotTracks(c(zoomucscTrack, splTrack, sTrack, alTrack, zoomcoverTrack), from = pos - 28, to = pos + 16, chromosome = chr)
##########

#
plotTracks(c(sTrack, splTrack, sTrack, alTrack, zoomcoverTrack), from = pos - 28, to = pos + 16, chromosome = chr)

##




##



start(ucscTrack)
transcript(ucscTrack)

range(ucscTrack)
#


# Data track
scores <- read.xlsx(scores_xlsx)
splTrack <- DataTrack(
  range = scores, chromosome = chr, genome = gen,name="∆ score", 
  cex.title = 1.2, background.title = "#F8ACAC", type = "histogram", 
  baseline = 0, ylim = c(-1, 1), lwd.baseline = 1,
  yTicksAt = c(-1.0, -0.5, 0, 0.5, 1.0), groups = c("AG", "AL", "DG", "DL"),
  col = c("#6088C6", "#EB8686", "#73D0C2", "#ED8D49"), legend = TRUE, 
  cex.legend = 1, size = 30
)

options(ucscChromosomeNames=FALSE) 


plotTracks(
  chromosome = chr,
  c(iTrack, ucscTrack, sTrack, alTrack, sashimiTrack), 
  # sashimiFilter = introns, sashimiFilterTolerance = 50,
  to = to + 300, from = from - 8000
)

##

introns <- GRanges("chrX", IRanges(
  start = c(154032556, 154092308, 154092308, 154032558, 154031451),
  end = c(154039791, 154097603, 154097619, 154092183, 154032206)))

introns <- GRanges("chrX", IRanges(
  start = c(154032556, 154092308, 154032558, 154031451),
  end = c(154039791, 154097603, 154092183, 154032206)))

introns <- GRanges("chrX", IRanges(
  start = c(154032556, 154092308, 154032558, 154032558),
  end = c(154039791, 154097603, 154097619, 154092183)))


#