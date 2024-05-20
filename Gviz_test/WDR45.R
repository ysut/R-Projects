BiocManager::install("txdbmaker")

library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(openxlsx)
library(rtracklayer)
library(GenomicRanges)
library(txdbmaker)

options(Gviz.ucscUrl="http://genome-asia.ucsc.edu/cgi-bin/")

session <- browserSession() # Open a browser session
genome(session) <- "hg19" # Set the genome to hg19

query <- ucscTableQuery(session, "NCBI RefSeq", 
                        GRangesForUCSCGenome(gen, chr, IRanges(from,to)))

query <- ucscTableQuery(session, "Conservation", 
                        GRangesForUCSCGenome(gen, chr, IRanges(from,to)))

tableNames(query)

ucscTrack <- UcscTrack(genome = gen, 
                       chromosome = chr, 
                       track = "NCBI RefSeq", 
                       from = from, 
                       to = to, 
                       trackType = "GeneRegionTrack", 
                       rstarts = "exonStarts", 
                       rends = "exonEnds", 
                       gene ="name", 
                       symbol = "name2", 
                       transcript = "name", 
                       strand = "strand" 
)


conservationTrack <- UcscTrack(genome = gen,
                               chromosome = chr,
                               range = range,
                               track = "Conservation",
                               table = "phyloP100wayAll",
                               trackType = "DataTrack",
                               start = "start", end = "end", 
                               data = "score",type = "hist", window = "auto",
                               col.histogram = "darkblue",fill.histogram = "darkblue",
                               ylim = c(-3.7, 4),name = "Conservation")


plotTracks(list(seqtrack, ucscTrack), from = pos-10, to = pos+10)



##### ==================================================================== #####

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
  

# Make GRanges object from gff3 file ===========================================

gff_file<-system.file("extdata", "GFF3_files", "fromucsc.gff3", package="txdbmaker")
txdb <- makeTxDbFromGFF(gff_file, format = "gff3")

txdb


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
  genome = gen, trackType = "GeneRegionTrack", fill = "#665990", name = "MECP2",
  track = "NCBI RefSeq", table = "ncbiRefSeq",
  rstarts = "exonStarts", rends = "exonEnds", gene = "name", symbol = "name2",  
  transcript = "name", strand = "strand"
  chromosome = mecp2_chr, from = mecp2_from, to = mecp2_to, 
  )


displayPars(refgenetrack) <- list(
  background.title = "#665990", cex = 1, transcriptAnnotation = "transcript",
  fontsize = 12)

transcript(refgenetrack) <- c("MANE", "MANE","MANE","MANE","MANE","MANE","MANE","MANE","MANE","MANE","MANE",
                              "MANE plus clinical", "MANE plus clinical","MANE plus clinical","MANE plus clinical","MANE plus clinical","MANE plus clinical",
                              "MANE plus clinical","MANE plus clinical","MANE plus clinical","MANE plus clinical","MANE plus clinical","MANE plus clinical")
transcript(refgenetrack)
initialize(refgenetrack)

track_list = list(itrack, axtrack, refgenetrack)
ht = HighlightTrack(
  trackList = track_list, alpha = 0.5, inBackground = FALSE, width = 0,
  chromosome = mecp2_chr, start = mecp2_pos
  
  )


plotTracks(list(ht, vartrack), from = pos - 20, to = pos + 20)
plotTracks(list(refgenetrack, vartrack), from = from - 5000, to = to + 20)

refgenetrack

# SpliceAI scores (DS_AG, DS_AL, DS_DG, DS_DL)
scores_xlsx <- '/Volumes/vol/work/Github/R-Projects/Gviz_test/splai_wdr45.xlsx'
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
seqtrack <- SequenceTrack(Hsapiens, chromosome = chr)

plotTracks(list(axtrack, seqtrack, vartrack, refgenetrack, splaitrack), 
           chromosome = chr, from = pos - 30, to = pos + 12)

