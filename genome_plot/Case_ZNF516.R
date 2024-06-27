library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(openxlsx)
library(rtracklayer)

## Tracks & Tables
session <- browserSession()
genome(session) <- "hg19"
query <- ucscTableQuery(session, "Fix Patches")
trackNames(query)
tableNames(query)

# Conncecting options 
# options(Gviz.ucscUrl="http://genome-asia.ucsc.edu/cgi-bin/")

# Genome assembly
gen <- "hg19"

znf516_chr  <- "chr18"
znf516_from <- 74069637
znf516_to   <- 74207198
znf516_pos  <- 74082552
znf516_bam  <- "bams/Sample_20287.recal.18_74069637-74207198.rh.bam"
znf516_xlsx <- "spliceai_scores/splai_znf516.xlsx"

chr  <- znf516_chr
from <- znf516_from
to   <- znf516_to
pos  <- znf516_pos
input_bam <- znf516_bam
scores_xlsx <- znf516_xlsx


## CCRs
ccrs_bw <- "bigwig/ccrs.autosomes.v2.20180420.bw"
grangesData <- import(ccrs_bw, format = "bigWig")
highScoreData <- grangesData[grangesData$score >= 95]
mediumScoreData <- grangesData[grangesData$score >= 90 & grangesData$score < 95]
lowScoreData <- grangesData[grangesData$score < 90]

highScoreTrack <- DataTrack(
  range = highScoreData, type = "histo", chromosome = chr, genome = gen,
  name = "CCRs", fill = "#EF8875", ylim = c(0, 100)
)

mediumScoreTrack <- DataTrack(
  range = mediumScoreData, type = "histo", chromosome = chr, genome = gen,
  name = "CCRs", fill = "#F7CDB6", ylim = c(0, 100), 
  fontsize = 10, size = 24, cex.title = 1, cex.axis = 0.5
)
lowScoreTrack <- DataTrack(
  range = lowScoreData, type = "histo", chromosome = chr, genome = gen,
  name = "CCRs (%tile)", fill = "lightgrey", ylim = c(0, 100), alpha = 0.5
)
ov = OverlayTrack(
  list(mediumScoreTrack, lowScoreTrack, highScoreTrack),
  name = "CCRs", background.title = "#909090"
  )

plotTracks(
  list(ov), from = from, to = to
)


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
  background.title = "#404040", fill ="#665990", cex = 1, cex.group = 1.2,
  transcriptAnnotation = "transcript",
  fontsize = 3, size = 0.2, cex.title = 2.5,
  rotation.title = 0
)


# Highlight the variant position
ht <- HighlightTrack(list(ucscTrack), alpha = 0.5, inBackground = FALSE, 
  chromosome = chr,
  start = c(pos - 80, 74074450), width = c(160, 100)
)

# Wide view with highlighted variant position
plotTracks(
  list(iTrack, ht, axTrack), from = from - 22000, to = to + 1000)

# ideogram
iTrack <- IdeogramTrack(
  genome = gen, chromosome = chr, 
  cex = 1.2, bevel = 1, showId = FALSE, size = 0.3, cex.bands = 1.2
)

### 2nd wide
# Highlight the variant position
ht2 <- HighlightTrack(list(ucscTrack, ov), alpha = 0.5, inBackground = FALSE, 
                     chromosome = chr,
                     start = c(pos - 80, 74074450), width = c(160, 100)
)

# Wide view with highlighted variant position
plotTracks(
  list(ht2), from = 74070000, to = 74070000 + 14000)


# chr18:74074452-74074513 and chr18:74082483-74082524
idrTrack <- AnnotationTrack(
  genome = gen, width = c(60, 41), name = 'IDR', fill = "darkgrey", 
  background.title = "#404040", shape = "box", group = "IDR", 
  just.group = "right", showOverplotting = TRUE, showTitle = TRUE,
  chromosome = chr, start = c(74074452, 74082483), size = 12, 
  rotation.title = 90, fontsize = 5, cex.title = 2.5
)

##########################
##### Zoomed in view #####
##########################
# axis
axTrack <- GenomeAxisTrack(
  add35 = FALSE, add53 = FALSE, exponent = 0, fontcolor = "#383838",
  fontsize = 12, size = 16, distFromAxis = 1.5
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
  name = "Gene", cex.title = 0.5, fontsize = 4, cex.axis = 1, rotation.title = 0,
  size = 6
)

##
fontcolor = c(A = "#49A190", C = "#0072B2", G = "#1E181B", T = "red", N = 1)
sTrack <- SequenceTrack(
  Hsapiens, fontcolor = fontcolor, fontsize = 6, 
  size = 5, 
  ) 

##
scores <- read.xlsx(scores_xlsx)
splTrack <- DataTrack(
  range = scores, chromosome = chr, genome = gen, name="∆ score", 
  cex.title = 1.0, background.title = "#F8ACAC", type = "histogram", 
  baseline = 0, col.baseline = "#838383", lty.baseline = 2, lwd.baseline = 2,
  ylim = c(-1, 1), lwd.baseline = 2,
  yTicksAt = c(-1.0, -0.5,0 , 0.5, 1.0), groups = c("AG", "AL", "DG", "DL"),
  col = c("#6088C6", "#EF8875", "#49A190", "#ED8D49"), legend = FALSE, 
  cex.legend = 2, cex.title = 1, fontsize = 12, cex.axis = 1, 
  size = 55
  )

##
options(ucscChromosomeNames=FALSE) 
alTrack <- AlignmentsTrack(
  input_bam, isPaired = TRUE, genome = "hg19", chromosome = chr, 
  background.title = "darkgrey", type = c("pileup"), 
  stacking = "squish", name = "WES reads", showAxis = FALSE,
  max.height = 10, min.height = 0.1, lwd.mismatch = 0.1, alpha.reads = 0.75, 
  col.axis = "darkgray", cex.title = 1, fontsize = 10, cex.axis = 1, 
  size = 25
)

plotTracks(c(axTrack, zoomucscTrack, splTrack, sTrack, alTrack, ov, idrTrack), 
           from = pos - 75, to = pos + 7, chromosome = chr
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
  col.axis = "darkgray", cex.title = 1, fontsize = 14, cex.axis = 1, 
  size = 20
)

sTrack <- SequenceTrack(
  Hsapiens, fontcolor = fontcolor, fontsize = 12, 
  size = 5
) 

ptcTrack <- AnnotationTrack(
    genome = gen, width = 2, name = 'PTC', fill = "red", 
    background.title = "#404040", shape = "box", group = "PTC", 
    just.group = "right", showOverplotting = TRUE, showTitle = FALSE,
    chromosome = chr, start = 74074496, size = 5, rotation.title = 0, fontsize = 8
)

plotTracks(c(axTrack, zoomucscTrack, ptcTrack, sTrack, alTrack, ov, idrTrack), 
           from = 74074486, to = 74074515, chromosome = chr
)