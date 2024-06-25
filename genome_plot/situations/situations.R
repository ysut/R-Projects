library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(openxlsx)

# Genome assembly and Chr.
gen <- "hg19"
chr <- "chr1"

###  Axis Track  ###
axTrack <- GenomeAxisTrack(
  add35 = FALSE, add53 = FALSE, exponent = 0, fontcolor = "#383838",
  fontsize = 10, labelPos = "above", size = 3, distFromAxis = 1.2
)

###  Dummy scores for overlap track ###
dummyscores <- read.xlsx("dummy_scores.xlsx")
dummyTrack <- DataTrack(
  range = dummyscores, chromosome = chr, genome = gen, name="SpliceAI ∆Score", 
  background.title = "#F8ACAC", type = "histogram", 
  baseline = 0, col.baseline = "#838383", lty.baseline = 1, lwd.baseline = 1,
  ylim = c(-1, 1),
  yTicksAt = c(-1.0, 0, 1.0), groups = c("AG", "AL", "DG", "DL"),
  col = c("#6088C6", "#EF8875", "#49A190", "#ED8D49"), legend = FALSE, 
  cex.title = 1.5, 
  fontsize = 12, cex.axis = 1.5, 
  size = 1
)

#####################
#1. Pseudo exon activation
#2. Exon skipping
#3. Whole intron retention
#4. Partial exon deletion
#5. Partial intron retention

#1. Pseudo exon activation
###  Gene model  ###
gene_model_xlsx <- "situations.xlsx"
genemodel <- read.xlsx(gene_model_xlsx)
gmtrack <- GeneRegionTrack(
  genemodel, name = "Gene Model", background.title = "#665990", fill = "#383838",
  chromosome = chr, genome = "hg19", size = 0.1, alpha = 0.70)

pseudoexon_xlsx <- "pseudoexon.xlsx"
pseudoexon_scores <- read.xlsx(pseudoexon_xlsx)
pseudoex_splaiTrack <- DataTrack(
  range = pseudoexon_scores, chromosome = chr, genome = gen, name="SpliceAI ∆Score", 
  background.title = "#F8ACAC", type = "histogram", 
  ylim = c(-1, 1), lwd.baseline = 2,
  yTicksAt = c(-1.0, 0, 1.0), groups = c("AG", "AL", "DG", "DL"),
  col = c("#6088C6", "#EF8875", "#49A190", "#ED8D49"), legend = FALSE, 
  cex.legend = 2, cex.title = 1.5, 
  fontsize = 12, cex.axis = 1.5, 
  size = 1
)


pseudo_xlsx <- "pseudo.xlsx"
pseudomodel <- read.xlsx(pseudo_xlsx)
pseudoexonTrack <- GeneRegionTrack(
  pseudomodel, name = "Gene Model", background.title = "#665990", fill = "red",
  chromosome = chr, genome = "hg19", size = 0.3, alpha = 0.3)




### SpliceAI
exonskip_xlsx <- "skip.xlsx"
exonskip_scores <- read.xlsx(exonskip_xlsx)
skipTrack <- DataTrack(
  range = exonskip_scores, chromosome = chr, genome = gen, name="SpliceAI ∆Score", 
  background.title = "#F8ACAC", type = "histogram", 
  baseline = 0, col.baseline = "#838383", lty.baseline = 2, lwd.baseline = 2,
  ylim = c(-1, 1), lwd.baseline = 2,
  yTicksAt = c(-1.0, 0, 1.0), groups = c("AG", "AL", "DG", "DL"),
  col = c("#6088C6", "#EF8875", "#49A190", "#ED8D49"), legend = FALSE, 
  cex.legend = 2, cex.title = 1.5, 
  fontsize = 9, cex.axis = 1.5, 
  size = 1.4
)

plotTracks(c(axTrack, gmtrack, skipTrack), from = 980, to = 1170)



ov <- OverlayTrack(
  list(dummyTrack, gmtrack, pseudoex_splaiTrack, pseudoexonTrack), background.title = "#665990"
)
plotTracks(c(axTrack, ov), from = 980, to = 1170)


##
ov <- OverlayTrack(
  list(gmtrack, pseudoex_splaiTrack), background.title = "#665990"
)
plotTracks(c(axTrack, ov, pseudoexonTrack), from = 980, to = 1170)

