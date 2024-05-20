library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(openxlsx)

# Genome assembly
gen <- "hg19"

# Gene physical location
chr <- "chrX"
from <- 48932101
to   <- 48958116

# Variant position
pos  <- 48933525

# Ideogram
itrack <- IdeogramTrack(
  genome = gen, chromosome = chr, cex = 0.8, bevel = 1, showId = FALSE
  )

# Axis track
axtrack <- GenomeAxisTrack(
  add35 = F, add53 = F, exponent = 0, fontcolor = "#383838", 
  fontsize = 16, lbaelPos = "below"
  )

# Fetch gene information from UCSC
refgenetrack <- UcscTrack(
  genome = gen, chromosome = chr, track = "NCBI RefSeq", table = "ncbiRefSeq",
  from = from, to = to, trackType = "GeneRegionTrack", rstarts = "exonStarts", 
  rends = "exonEnds", gene = "name", symbol = "name",  transcript = "name", 
  strand = "strand", fill = "#665990", name = "WDR45", geneSymbols = TRUE
  )

displayPars(refgenetrack) <- list(
  background.title = "#665990", cex = 2, transcriptAnnotation = "transcript",
  fontsize = 15)

track_list = list(itrack, axtrack, refgenetrack)
ht = HighlightTrack(trackList = track_list, chromosome = chr, 
                    start = pos, width = 0,
                    alpha = 0.5, inBackground = FALSE
                    )

plotTracks(list(ht, vartrack), from = pos - 20, to = pos + 20)

refgenetrack

# SpliceAI scores (DS_AG, DS_AL, DS_DG, DS_DL)
scores_xlsx <- '/Volumes/vol/work/Github/R-Projects/Gviz_test/wdr45.xlsx'
scores <- read.xlsx(scores_xlsx)
splaitrack <- DataTrack(
  range = scores, chromosome = chr, genome="hg19", name="SpliceAI delta score", 
  cex.title = 1.0, background.title = "#F8ACAC", type = "histogram", 
  baseline = 0, ylim = c(-1, 1), lwd.baseline = 1, 
  yTicksAt = c(-1.0, -0.2, 0, 0.2, 1.0), groups = c("AG", "AL", "DG", "DL"),
  col = c("#6088C6", "#EB8686", "#73D0C2", "#ED8D49"),
  cex.legend = 0.7, box.legend = FALSE)

# Variant position
vartrack <- AnnotationTrack(
  genome = gen, chromosome = chr, start = pos, width = 0, name = 'Variant Pos.',
  fill = "#383838", background.title = "#383838", shape = "box", group = "Variant",
  just.group = "right", showOverplotting = TRUE)

# Sequence track
seqtrack <- SequenceTrack(Hsapiens)

plotTracks(list(axtrack, seqtrack, vartrack, refgenetrack, splaitrack), 
           chromosome = chr, from = pos - 80, to = pos + 80)

