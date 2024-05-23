library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(openxlsx)
library(txdbmaker)

# Load the data

# Genome assembly
gen <- "hg19"

# Gene physical location
mecp2_chr  <- "chrX"
mecp2_from <- 48932101
mecp2_to   <- 48958116

col2a1_chr  <- "chr12"
col2a1_from <- 48366750
col2a1_to   <- 48398259

pdha1_chr   <- "chrX"
pdha1_from  <- 19362045
pdha1_to    <- 19379836

jakmip1_chr  <- "chr4"
jakmip1_from <- 6027926
jakmip1_to   <- 6202276

ubn1_chr   <- "chr16"
ubn1_from  <- 4897482
ubn1_to    <- 4932402

znf516_chr  <- "chr18"
znf516_from <- 74069637
znf516_to   <- 74207198

nfe2l1_chr  <- "chr17"
nfe2l1_from <- 46125721
nfe2l1_to   <- 46138907

# Variant position
mecp2_pos   <- 48933525
col2a1_pos  <- 48387611
peha1_pos   <- 19373601
jakmip1_pos <- 6086572
ubn1_pos    <- 4918905
znf516_pos  <- 74082552
nfe2l1_pos  <- 46134392

# Gene names
options(Gviz.ucscUrl="http://genome-asia.ucsc.edu/cgi-bin/")
# Open a UCSC session
session <- browserSession("UCSC")
genome(session) <- "hg19"
# List available tracks
tracks <- trackNames(session)
print(tracks)

## Variant info
chr  <- col2a1_chr
from <- col2a1_from
to   <- col2a1_to
pos  <- col2a1_pos
scores_xlsx <- 'splai_col2a1.xlsx'

chr  <- mecp2_chr
from <- mecp2_from
to   <- mecp2_to
pos  <- mecp2_pos
scores_xlsx <- 'splai_mecp2.xlsx'

##################################################
# For wide view
##################################################
# ideogram
iTrack <- IdeogramTrack(
  genome = gen, chromosome = chr, cex = 2, bevel = 1, showId = TRUE
  )

# axis
axTrack <- GenomeAxisTrack(
  add35 = FALSE, add53 = FALSE, exponent = 0, fontcolor = "#383838",
  fontsize = 16, labelPos = "above"
  )

ucscTrack <- UcscTrack(
  genome = gen, track = "NCBI RefSeq", trackType = "GeneRegionTrack",
  rstarts = "exonStarts", rends = "exonEnds", gene ="name", symbol = "name2", 
  transcript = "name", strand = "strand", table = "ncbiRefSeqSelect",
  chromosome = chr, from = from, to = to, name = "COL2A1"
)

displayPars(ucscTrack) <- list(
  background.title = "#665990", fill ="#665990", cex = 1, 
  transcriptAnnotation = "transcript", fontsize = 12
  )

for (i in 1:length(start(ucscTrack))) {
  start(ucscTrack)[i] <- start(ucscTrack)[i] + 1
}

# Highlight the variant position
ht <- HighlightTrack(ucscTrack, alpha = 0.5, inBackground = FALSE, 
  chromosome = chr, start = pos -100, width = 200
)

# Wide view with highlighted variant position
plotTracks(
  list(iTrack, axTrack, ht), from = from - 1500, to = to + 500, type = "none"
  )


##################################################
# For zoomed in view
##################################################

# Sequence track
sTrack <- SequenceTrack(Hsapiens, chromosome = chr) 

# Variant position
varTrack <- AnnotationTrack(
  genome = gen, width = 0, name = 'Var. Pos.', fill = "#383838", 
  background.title = "#383838", shape = "box", group = "Variant",
  just.group = "right", showOverplotting = TRUE, 
  chromosome = chr, start = pos, 
  )

# Data track
scores <- read.xlsx(scores_xlsx)
splTrack <- DataTrack(
  range = scores, chromosome = chr, genome="hg19", name="SpliceAI ∆", 
  cex.title = 1.5, background.title = "#F8ACAC", type = "histogram", 
  baseline = 0, ylim = c(-1, 1), lwd.baseline = 1, 
  yTicksAt = c(-1.0, -0.5, 0, 0.5, 1.0), groups = c("AG", "AL", "DG", "DL"),
  col = c("#6088C6", "#EB8686", "#73D0C2", "#ED8D49"),
  cex.legend = 0.7, box.legend = FALSE)


# Zoom in to the variant position
plotTracks(list(sTrack, varTrack, ucscTrack, splTrack), chromosome = chr, from = pos - 50, to = pos + 50)

# Conservation track
# STEP 1: Fetach conservation data
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


# Zoom in to the variant position
plotTracks(list(sTrack, varTrack, ucscTrack, splTrack, absplice, phylop100, gerp), from = pos - 32, to = pos + 45)
plotTracks(list(axTrack, sTrack, varTrack, ucscTrack, splTrack, absplice, consot), from = pos - 50, to = pos + 50)
