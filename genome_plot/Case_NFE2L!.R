library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(openxlsx)

# Conncecting options
options(Gviz.ucscUrl="http://genome-asia.ucsc.edu/cgi-bin/")

# Genome assembly
gen <- "hg19"

nfe2l1_chr  <- "chr17"
nfe2l1_from <- 46125721
nfe2l1_to   <- 46138907
nfe2l1_pos  <- 46134392
nfe2l1_bam  <- "bams/Sample_12988.recal.17_46125721-46138907.rh.bam"
nfe2l1_xlsx <- "spliceai_scores/splai_nfe2l1.xlsx"

chr  <- nfe2l1_chr
from <- nfe2l1_from
to   <- nfe2l1_to
pos  <- nfe2l1_pos
input_bam <- nfe2l1_bam
scores_xlsx <- nfe2l1_xlsx

## CCRs
ccrs_bw <- "bigwig/ccrs.autosomes.v2.20180420.bw"
grangesData <- import(ccrs_bw, format = "bigWig")
highScoreData <- grangesData[grangesData$score >= 95]
mediumScoreData <- grangesData[grangesData$score >= 90 & grangesData$score < 95]
lowScoreData <- grangesData[grangesData$score < 90]

highScoreTrack <- DataTrack(
  range = highScoreData, type = "histo", chromosome = chr, genome = gen,
  name = "CCR", fill = "#EF8875", ylim = c(0, 100)
)
mediumScoreTrack <- DataTrack(
  range = mediumScoreData, type = "histo", chromosome = chr, genome = gen,
  name = "CCR", fill = "#F7CDB6", ylim = c(0, 100), 
  fontsize = 4, cex.title = 3, cex.asix = 3, size = 30,
  yTicksAt = c(0, 50 , 100)
)
lowScoreTrack <- DataTrack(
  range = lowScoreData, type = "histo", chromosome = chr, genome = gen,
  name = "CCR", fill = "#303030", ylim = c(0, 100), alpha = 0.3
)
ov = OverlayTrack(
  list(mediumScoreTrack, lowScoreTrack, highScoreTrack),
  name = "CCRs (%tile)", background.title = "#909090"
)
plotTracks(ov,  from = pos - 5, to = pos + 98, chromosome = chr)

#####################
### For wide view ###
#####################
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
  background.title = "#404040", fill ="#665990", cex = 1.1, cex.group = 1.2,
  transcriptAnnotation = "transcript",
  fontsize = 2.0, size = 0.6, cex.title = 4,
  rotation.title = 0
)

# Highlight the variant position
ht <- HighlightTrack(ucscTrack, alpha = 0.5, inBackground = FALSE, 
  chromosome = chr, 
  start = c(pos - 40), width = c(180)
)
# ideogram
iTrack <- IdeogramTrack(
  genome = gen, chromosome = chr, 
  cex = 1.2, bevel = 1, showId = FALSE, size = 0.5, cex.bands = 1.2
)
# Wide view with highlighted variant position
plotTracks(
  list(iTrack, ht, axTrack), from = from - 3000, to = to + 1000
)

##########################
##### Zoomed in view #####
##########################
# axis
axTrack <- GenomeAxisTrack(
  add35 = FALSE, add53 = FALSE, exponent = 0, fontcolor = "#383838",
  fontsize = 14, size = 16, distFromAxis = 1.5
)

zoomucscTrack <- UcscTrack(
  genome = gen, track = "NCBI RefSeq", trackType = "GeneRegionTrack",
  rstarts = "exonStarts", rends = "exonEnds", gene ="name", symbol = "name2", 
  transcript = "name", strand = "strand", table = "ncbiRefSeqSelect",
  chromosome = chr, from = pos - 1000, to = pos + 1000, name = "Gene"
)
for (i in 1:length(start(zoomucscTrack))) {
  start(zoomucscTrack)[i] <- start(zoomucscTrack)[i] + 1
}

displayPars(zoomucscTrack) <- list(
  background.title = "#404040", fill ="#665990", cex = 1.0, showTitle = FALSE,
  name = "Gene", cex.title = 3, fontsize = 4, cex.axis = 4, rotation.title = 0,
  size = 7
)

##
fontcolor = c(A = "#49A190", C = "#0072B2", G = "#1E181B", T = "red", N = 1)
sTrack <- SequenceTrack(
  Hsapiens, fontcolor = fontcolor, fontsize = 12, 
  size = 6, 
  ) 

##
scores <- read.xlsx(scores_xlsx)
splTrack <- DataTrack(
  range = scores, chromosome = chr, genome = gen, name="∆ score", 
  background.title = "#F8ACAC", type = "histogram", 
  baseline = 0, col.baseline = "#838383", lty.baseline = 2, lwd.baseline = 2,
  ylim = c(-1, 1), lwd.baseline = 2,
  yTicksAt = c(-1.0, -0.5,0 , 0.5, 1.0), groups = c("AG", "AL", "DG", "DL"),
  col = c("#6088C6", "#EF8875", "#49A190", "#ED8D49"), legend = FALSE, 
  cex.legend = 2, cex.title = 3, fontsize = 6, cex.axis = 2, 
  size = 55
  )

##
options(ucscChromosomeNames=FALSE) 
alTrack <- AlignmentsTrack(
  input_bam, isPaired = TRUE, genome = "hg19", chromosome = chr, 
  background.title = "darkgrey", type = c("pileup"), 
  stacking = "squish", name = "WES reads", showAxis = FALSE,
  max.height = 10, min.height = 1, lwd.mismatch = 0.1, alpha.reads = 0.75, 
  col.axis = "darkgray", cex.title = 3, fontsize = 6,
  size = 50
)
hl2 <- HighlightTrack(
  list(ptcTrack, sTrack, alTrack), 
  alpha = 0.5, inBackground = FALSE, chromosome = chr, 
  start = c(pos+46), width = c(8)
)
plotTracks(c(axTrack, zoomucscTrack, splTrack, hl2, ov), 
           from = pos - 8, to = pos + 130, chromosome = chr
)

s2Track <- SequenceTrack(
  Hsapiens, fontcolor = fontcolor, fontsize = 22, 
  size = 8, 
) 
plotTracks(c(zoomucscTrack, splTrack, ptcTrack, s2Track, alTrack), 
           from = pos + 46, to = pos + 55, chromosome = chr
)

#################
### PTC Track ###
#################

options(ucscChromosomeNames=FALSE)
alTrack <- AlignmentsTrack(
  input_bam, isPaired = TRUE, genome = "hg19", chromosome = chr, 
  background.title = "darkgrey", type = c("pileup"), 
  stacking = "squish", name = "WES reads", showAxis = FALSE,
  max.height = 10, min.height = 0.1, lwd.mismatch = 0.1, alpha.reads = 0.75, 
  col.axis = "darkgray", cex.title = 12, fontsize = 2, cex.axis = 3, 
  size = 20
)

sTrack <- SequenceTrack(
  Hsapiens, fontcolor = fontcolor, fontsize = 12, 
  size = 5
) 

ptcTrack <- AnnotationTrack(
    genome = gen, width = 2, name = 'PTC', fill = "red", 
    background.title = "#404040", shape = "box", group = "PTC", 
    just.group = "right", showOverplotting = TRUE,
    chromosome = chr, start = 46134441, size = 7, fontsize = 24,
    rotation.title = 0
)


plotTracks(c(axTrack, zoomucscTrack, ptcTrack, sTrack, alTrack), 
           from = 4920325, to = 4920350, chromosome = chr)
