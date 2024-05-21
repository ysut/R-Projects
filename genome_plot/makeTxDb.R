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


# Create the TxDb object

gffFile <- system.file("extdata", "GFF3_files", "col2a1.gff", package="GenomicFeatures")
txdb <- makeTxDbFromGFF(file=gffFile,
                        dataSource="partial gtf file for Tomatoes for testing",
                        organism="Solanum lycopersicum")

## TESTING GTF, this time specifying the chrominfo


chr  <- col2a1_chr
from <- col2a1_from
to   <- col2a1_to
pos  <- col2a1_pos

iTrack <- IdeogramTrack(
  genome = gen, chromosome = chr, cex = 2, bevel = 1, showId = TRUE
  )

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
  list(iTrack, axTrack, ht), from = from - 2000, to = to + 500, type = "none"
  )


##################################################
# For zoomed in view
##################################################

# Sequence track
sTrack <- SequenceTrack(Hsapiens, chromosome = chr) 

# Variant position
varTrack <- AnnotationTrack(
  genome = gen, width = 0, name = 'Variant Pos.', fill = "#383838", 
  background.title = "#383838", shape = "box", group = "Variant",
  just.group = "right", showOverplotting = TRUE, 
  chromosome = chr, start = pos, 
  )

# Data track



# Zoom in to the variant position
plotTracks(list(varTrack, sTrack, ucscTrack), from = pos - 50, to = pos + 50, type = "none")




