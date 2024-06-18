library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(openxlsx)
library(rtracklayer)

# Conncecting options
options(Gviz.ucscUrl="http://genome-asia.ucsc.edu/cgi-bin/")

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
  name = "CCRs (%tile)", fill = "#EF8875", ylim = c(0, 100)
)

mediumScoreTrack <- DataTrack(
  range = mediumScoreData, type = "histo", chromosome = chr, genome = gen,
  name = "CCRs (%tile)", fill = "#F7CDB6", ylim = c(0, 100), 
  fontsize = 1.2, size = 0.5, cex.title = 10,
)

lowScoreTrack <- DataTrack(
  range = lowScoreData, type = "histo", chromosome = chr, genome = gen,
  name = "CCRs (%tile)", fill = "lightgrey", ylim = c(0, 100), alpha = 0.3
)
ov = OverlayTrack(
  list(mediumScoreTrack, lowScoreTrack, highScoreTrack),
  name = "CCRs (%tile)", background.title = "#909090"
  )

plotTracks(
  list(iTrack, ht, ov, axTrack), from = 74153000, to = 74156000
)

#####################
# Disordered region
disordered_bw <- "bigwig/disordered_regions.bw"


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
  fontsize = 1.2, size = 0.5, cex.title = 10,
  rotation.title = 0
)

# Highlight the variant position
ht <- HighlightTrack(ucscTrack, alpha = 0.5, inBackground = FALSE, 
  chromosome = chr, 
  start = c(pos - 80, 4920300), width = c(160, 100)
)
# ideogram
iTrack <- IdeogramTrack(
  genome = gen, chromosome = chr, 
  cex = 1.2, bevel = 1, showId = FALSE, size = 0.3, cex.bands = 1.2
)
# Wide view with highlighted variant position
plotTracks(
  list(iTrack, ht, ov, axTrack), from = from - 15000, to = to + 1000
)


##########################
##### Zoomed in view #####
##########################
# axis
axTrack <- GenomeAxisTrack(
  add35 = FALSE, add53 = FALSE, exponent = 0, fontcolor = "#383838",
  fontsize = 20, size = 12, distFromAxis = 1.5
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
  name = "Gene", cex.title = 6, fontsize = 4, cex.axis = 4, rotation.title = 0,
  size = 6
)

##
fontcolor = c(A = "#49A190", C = "#0072B2", G = "#1E181B", T = "red", N = 1)
sTrack <- SequenceTrack(
  Hsapiens, fontcolor = fontcolor, fontsize = 12, 
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
  cex.legend = 2, cex.title = 6, fontsize = 6, cex.axis = 4, 
  size = 55
  )


##
options(ucscChromosomeNames=FALSE) 
alTrack <- AlignmentsTrack(
  input_bam, isPaired = TRUE, genome = "hg19", chromosome = chr, 
  background.title = "darkgrey", type = c("pileup"), 
  stacking = "squish", name = "WES reads", showAxis = FALSE,
  max.height = 10, min.height = 0.1, lwd.mismatch = 0.1, alpha.reads = 0.75, 
  col.axis = "darkgray", cex.title = 6, fontsize = 6, cex.axis = 4, 
  size = 25
)

plotTracks(c(axTrack, zoomucscTrack, splTrack, sTrack, alTrack), 
           from = pos - 100, to = pos + 5, chromosome = chr
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
    chromosome = chr, start = 4920336, size = 5, rotation.title = 0, fontsize = 24
)

plotTracks(c(axTrack, zoomucscTrack, ptcTrack, sTrack, alTrack), 
           from = 4920325, to = 4920350, chromosome = chr
)


##### CCRs


plotTracks(c(sTrack, ccrTrack), from = pos -100, to = pos +100, chromosome = chr)

