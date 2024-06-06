library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(openxlsx)
library(txdbmaker)

# Genome assembly
gen <- "hg19"

# Case 4
jakmip1_chr  <- "chr4"
jakmip1_from <- 6027926
jakmip1_to   <- 6202276
jakmip1_pos <- 6086572
jakmip1_bam <- "bams/Sample_21599.recal.4_6027926-6202276.rh.bam"
jakmip1_xlsx <- "splai_jakmip1.xlsx"


# Gene names
options(Gviz.ucscUrl="http://genome-asia.ucsc.edu/cgi-bin/")

# Case 4
chr  <- jakmip1_chr
from <- jakmip1_from
to   <- jakmip1_to
pos  <- jakmip1_pos
scores_xlsx <- jakmip1_xlsx
input_bam <- jakmip1_bam

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
  transcript = "name", strand = "strand", table = "ncbiRefSeqSelect",
  chromosome = chr, from = from, to = to, name = "Gene"
)

for (i in 1:length(start(ucscTrack))) {
  start(ucscTrack)[i] <- start(ucscTrack)[i] + 1
}

displayPars(ucscTrack) <- list(
  background.title = "#665990", fill ="#665990", cex = 4, 
  transcriptAnnotation = "transcript", fontsize = 12, size = 0.5,
  cex.group = 0.8, showId = TRUE
)

# Highlight the variant position
ht <- HighlightTrack(ucscTrack, alpha = 0.5, inBackground = FALSE, 
  chromosome = chr, start = pos - 100, width = 200
)

# Wide view with highlighted variant position
plotTracks(
  list(iTrack, ht, axTrack), from = from - 12000, to = to + 2000, type = "none"
)

##################################################
# For zoomed in view
##################################################
# Sequence track
sTrack <- SequenceTrack(Hsapiens) 

# Axis 
axTrack <- GenomeAxisTrack(
  add35 = FALSE, add53 = FALSE, exponent = 0, fontcolor = "#383838",
  fontsize = 16, labelPos = "above", size = 8, distFromAxis = 1.2
)
plotTracks(list(axTrack), from = pos - 10, to = pos + 10)

start(axTrack)
#


# Variant position
varTrack <- AnnotationTrack(
  genome = gen, width = 0, name = 'Var.', fill = "#383838", 
  background.title = "#383838", shape = "box", group = "Variant",
  just.group = "right", showOverplotting = TRUE, cex.title = 1.2,
  chromosome = chr, start = pos, size = 3
)

ucscTrack <- UcscTrack(
  genome = gen, track = "NCBI RefSeq", trackType = "GeneRegionTrack",
  rstarts = "exonStarts", rends = "exonEnds", gene ="name", symbol = "name2", 
  transcript = "name", strand = "strand", table = "ncbiRefSeqSelect",
  chromosome = chr, from = pos - 200, to = pos + 100, name = "",
  background.title = "#665990", fill ="#665990", cex.title = 1.0,
  exonAnnotations = "transcript", fontsize = 12, size = 3
)
for (i in 1:length(start(ucscTrack))) {
  start(ucscTrack)[i] <- start(ucscTrack)[i] + 1
}

# Data track
scores <- read.xlsx(scores_xlsx)
splTrack <- DataTrack(
  range = scores, chromosome = chr, genome="hg19", name="∆ score", 
  cex.title = 1.2, background.title = "#F8ACAC", type = "histogram", 
  baseline = 0, ylim = c(-1, 1), lwd.baseline = 1,
  yTicksAt = c(-1.0, -0.5, 0, 0.5, 1.0), groups = c("AG", "AL", "DG", "DL"),
  col = c("#6088C6", "#EB8686", "#73D0C2", "#ED8D49"), legend = TRUE, 
  cex.legend = 1, size = 30
)

options(ucscChromosomeNames=FALSE) 

alTrack <- AlignmentsTrack(
  input_bam, isPaired = TRUE, genome = "hg19", chromosome = chr, 
  background.title = "darkgrey", type = c("pileup"), 
  stacking = "squish", minCoverageHeight = 20, cex.title = 1.2,
  name = "BAM data", coverageHeight = 0.05, showAxis = FALSE,
  max.height = 8, min.height = 1, lwd.mismatch = 0.1, alpha.reads = 0.75, size = 16,
  col.axis = "darkgray"
)

plotTracks(
  list(axTrack, ucscTrack, splTrack, sTrack, alTrack),
  # list(axTrack, sTrack, varTrack, ucscTrack, splTrack, alTrack),
  chromosome = chr, from = pos-1180, to = pos+40
)


plotTracks(
  list(axTrack, ucscTrack, splTrack, sTrack, alTrack),
  # list(axTrack, sTrack, varTrack, ucscTrack, splTrack, alTrack),
  chromosome = chr, from = pos-1165, to = pos-1145
)


plotTracks(
  list(axTrack, ucscTrack, splTrack, sTrack, alTrack),
  # list(axTrack, sTrack, varTrack, ucscTrack, splTrack, alTrack),
  chromosome = chr, from = pos-10, to = pos + 10
)


mdTrack <- DataTrack(
  range = GRanges(
    seqnames = rep(chroms, c(10, 40, 20, 100)),
                  ranges = IRanges(start = c(seq(1, 100, len = 10),
                                             seq(1, 400, len = 40), 
                                             seq(1, 200, len = 20),
                                             seq(1, 1000, len = 100)), 
                                   width = 9), values = runif(170)),
  data = "values", chromosome = "chr1", genome = "mm9", name = "bar")


chroms <- c("chr1", "chr2", "chr3", "chr4")
range = GRanges(seqnames = rep(chroms, c(10, 40, 20, 100)),
                ranges = IRanges(start = c(seq(1, 100, len = 10),
                                           seq(1, 400, len = 40), 
                                           seq(1, 200, len = 20),
                                           seq(1, 1000, len = 100)), 
                                 width = 9), values = runif(170))

range


splTrack <- DataTrack(
  range = scores, chromosome = chr, genome="hg19", name="∆ score", 
  cex.title = 1.2, background.title = "#F8ACAC", type = "histogram", 
  baseline = 0, ylim = c(-1, 1), lwd.baseline = 1,
  yTicksAt = c(-1.0, -0.5, 0, 0.5, 1.0), groups = c("AG", "AL", "DG", "DL"),
  col = c("#6088C6", "#EB8686", "#73D0C2", "#ED8D49"), legend = TRUE, 
  cex.legend = 1, size = 30
)







##
######################## Conservation track ####################################
# Open a UCSC session
session <- browserSession("UCSC")
genome(session) <- "hg19"
# List available tracks
tracks <- trackNames(session)
print(tracks)

query1 <- ucscTableQuery(session, track = "Conservation", table = "phastCons100way")
query2 <- ucscTableQuery(session, track = "Conservation", table = "phyloP100wayAll")
query3 <- ucscTableQuery(session, track = "AbSplice Scores", table = "abSplice")
query4 <- ucscTableQuery(session, track = "GERP", table = "allHg19RS_BW")

range(query1) <- GRanges(seqnames = chr, ranges = IRanges(start = from, end = to))
range(query2) <- GRanges(seqnames = chr, ranges = IRanges(start = from, end = to))
range(query3) <- GRanges(seqnames = chr, ranges = IRanges(start = pos - 200, end = pos + 200))
range(query4) <- GRanges(seqnames = chr, ranges = IRanges(start = pos - 500, end = pos + 500))
data1 <- getTable(query1)
data2 <- getTable(query2)
data3 <- getTable(query3)
data4 <- getTable(query4)

# Check data only head
print(head(data1))
print(head(data2))
print(head(data3))
print(head(data4))

# STEP 2: Create a conservation track
data1$start <- as.numeric(data1$start)
data1$end <- as.numeric(data1$end)
data1$value <- as.numeric(data1$value)

data2$start <- as.numeric(data2$start)
data2$end <- as.numeric(data2$end)
data2$value <- as.numeric(data2$value)

data3$chromStart <- as.numeric(data3$chromStart)
data3$chromEnd <- as.numeric(data3$chromEnd)
data3$spliceABscore <- as.numeric(data3$spliceABscore)

data4$start <- as.numeric(data4$start)
data4$end <- as.numeric(data4$end)
data4$value <- as.numeric(data4$value)

# Reduce 1 bp from the end position in data4
data2$end <- data2$end - 1
data3$chromStart <- data3$chromStart + 1
data4$end <- data4$end - 1

# GRange objects
gr1 <- GRanges(
  seqnames = Rle(chr),
  ranges = IRanges(start = data1$start, end = data1$end),
  score = data1$value
)

gr2 <- GRanges(
  seqnames = Rle(chr),
  ranges = IRanges(start = data2$start, end = data2$end),
  score = data2$value
)

gr3 <- GRanges(
  seqnames = Rle(chr), 
  ranges = IRanges(start = data3$chromStart, end = data3$chromEnd),
  score = data3$spliceABscore
)

gr4 <- GRanges(
  seqnames = Rle(chr),
  ranges = IRanges(start = data4$start, end = data4$end),
  score = data4$value
)

phastcons100 <- DataTrack(
  range = gr1, genome = gen, chromosome = chr, data = data1$value,
  type = "hist", col.histogram = "darkblue", fill.histogram = "darkblue",
  ylim = c(0, 1), yTicksAt = c(0, 1.0), name = "phastsons100"
)

absplice <- DataTrack(
  range = gr3, genome = gen, chromosome = chr, data = "score",
  # col.histogram = "darkred", fill.histogram = "darkred",
  ylim = c(0, 0.5), yTicksAt = c(0, 0.1, 0.5), name = "AbSplice", type = "heatmap"
)
#
plotTracks(absplice, from = pos - 50, to = pos + 50)
#

phylop100 <- DataTrack(
  range = gr2, genome = gen, chromosome = chr, data = data2$value,
  type = "hist", col.histogram = "darkgreen", fill.histogram = "darkgreen",
  ylim = c(-8, 8), name = "PhyloP", alpha = 0.8
)

gerp <- DataTrack(
  range = gr4, genome = gen, chromosome = chr, data = data4$value,
  type = "hist", col.histogram = "darkred", fill.histogram = "darkred",
  ylim = c(-8, 8), name = "GERP", alpha = 0.8
)

consot <- OverlayTrack(
  list(phylop100, gerp), col = c("darkgreen", "darkred")
)