library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(openxlsx)

# Genome assembly and Chr.
gen <- "hg19"
chr <- "chr1"

###  Axis Track  ###
axTrack <- GenomeAxisTrack(
  add35 = FALSE, add53 = FALSE, exponent = 0, fontcolor = "#383838",
  fontsize = 10, labelPos = "above", size = 10, distFromAxis = 1.2
)

###  Dummy scores for overlap track ###
dummyscores <- read.xlsx("dummy_scores.xlsx")
dummyTrack <- DataTrack(
  range = dummyscores, chromosome = chr, genome = gen, name="SpliceAI\n∆Score", 
  background.title = "#F8ACAC", type = "histogram", ylim = c(-1.2, 1.2),
  yTicksAt = c(-1.0, 0, 1.0),
  groups = c("AG", "AL", "DG", "DL"),
  col = c("#4D4298", "#E60011", "#06A384", "#AB961D"), legend = FALSE, 
  cex.title = 1.0, 
  fontsize = 10, cex.axis = 1.0, 
  size = 2
)

#####################
#1. Pseudo exon activation
#2. Exon skipping
#3. Whole intron retention
#4. Partial exon deletion
#5. Partial intron retention

#1. Pseudo exon activation
gene_model_xlsx <- "01_pseudoexon/model_base.xlsx"
genemodel <- read.xlsx(gene_model_xlsx)
gmtrack <- GeneRegionTrack(
  genemodel, name = "Gene Model", background.title = "#665990", fill = "#383838",
  chromosome = chr, genome = "hg19", size = 0.1, alpha = 0.70)

pseudo_xlsx <- "01_pseudoexon/model_pseudoexon.xlsx"
pseudomodel <- read.xlsx(pseudo_xlsx)
pseudoexonTrack <- GeneRegionTrack(
  pseudomodel, name = "Gene Model", background.title = "#665990", fill = "red",
  chromosome = chr, genome = "hg19", size = 0.3, alpha = 0.2)

pseudoexon_xlsx <- "01_pseudoexon/splai_pseudoexon.xlsx"
pseudoexon_scores <- read.xlsx(pseudoexon_xlsx)
pseudoex_splaiTrack <- DataTrack(
  range = pseudoexon_scores, chromosome = chr, genome = gen, 
  background.title = "#F8ACAC", type = "histogram", 
  ylim = c(-1, 1), yTicksAt = c(-1, 0, 1),
  baseline = 0, col.baseline = "#838383", lty.baseline = 1, lwd.baseline = 1,
  groups = c("AG", "AL", "DG", "DL", "Variant"),
  col = c("#4D4298", "#E60011", "#06A384", "#AB961D", "black"), legend = FALSE, 
  size = 1
)

ov <- OverlayTrack(
  list(dummyTrack, gmtrack, pseudoex_splaiTrack, pseudoexonTrack), 
  background.title = "#838383"
)
plotTracks(c(axTrack, ov), from = 985, to = 1155)


#2. Exon skipping
genemodel <- read.xlsx("02_exskip/model_base_skip.xlsx")
gmtrack <- GeneRegionTrack(
  genemodel, fill = "#383838",
  chromosome = chr, genome = "hg19", size = 0.1, alpha = 0.70)

skipex <- read.xlsx("02_exskip/model_exonskip.xlsx")
skipexTrack <- GeneRegionTrack(
  skipex, fill = "#383838", chromosome = chr, genome = "hg19", 
  size = 0.3, alpha = 0.2)

exonskip_scores <- read.xlsx("02_exskip/splai_exonskip.xlsx")
exonskip_splaiTrack <- DataTrack(
  range = exonskip_scores, chromosome = chr, genome = gen, name="SpliceAI ∆Score", 
  background.title = "#F8ACAC", type = "histogram", 
  ylim = c(-1, 1), yTicksAt = c(-1, 0, 1),
  baseline = 0, col.baseline = "#838383", lty.baseline = 1, lwd.baseline = 1,
  groups = c("AG", "AL", "DG", "DL", "Var"),,
  col = c("#4D4298", "#E60011", "#06A384", "#AB961D", "black"), legend = FALSE, 
  size = 1
)

ov <- OverlayTrack(
  list(dummyTrack, gmtrack, exonskip_splaiTrack, skipexTrack), 
  background.title = "#838383"
)
plotTracks(c(axTrack, ov), from = 985, to = 1155)


#3. Whole intron retention
genemodel <- read.xlsx("03_intret/model_base.xlsx")
gmtrack <- GeneRegionTrack(
  genemodel, fill = "#383838",
  chromosome = chr, genome = "hg19", size = 0.1, alpha = 0.60)

intret <- read.xlsx("03_intret/model_intret.xlsx")
intretTrack <- GeneRegionTrack(
  intret, fill = "red",
  chromosome = chr, genome = "hg19", size = 0.3, alpha = 0.2)

intret_scores <- read.xlsx("03_intret/splai_whole_intret.xlsx")
intret_splaiTrack <- DataTrack(
  range = intret_scores, chromosome = chr, genome = gen, name="SpliceAI ∆Score", 
  background.title = "#F8ACAC", type = "histogram", 
  ylim = c(-1, 1), yTicksAt = c(-1, 0, 1),
  baseline = 0, col.baseline = "#838383", lty.baseline = 1, lwd.baseline = 1.0,
  groups = c("AG", "AL", "DG", "DL", "Var"),
  col = c("#4D4298", "#E60011", "#06A384", "#AB961D", "black"), legend = FALSE, 
  size = 1
)

ov <- OverlayTrack(
  list(dummyTrack, gmtrack, intret_splaiTrack, intretTrack), 
  background.title = "#838383"
)
plotTracks(c(axTrack, ov), from = 985, to = 1155)

#4. Partial exon deletion
genemodel <- read.xlsx("04_partexdel/model_base_exdel.xlsx")
genemodel <- read.xlsx("04_partexdel/model_base_exdel_2.xlsx")
gmtrack <- GeneRegionTrack(
  genemodel, fill = "#383838",
  chromosome = chr, genome = "hg19", size = 0.1, alpha = 0.60)

exdel <- read.xlsx("04_partexdel/model_exdel.xlsx")
exdel <- read.xlsx("04_partexdel/model_exdel_2.xlsx")
exdelTrack <- GeneRegionTrack(
  exdel, fill = "#383838",
  chromosome = chr, genome = "hg19", size = 0.3, alpha = 0.2)

exdel_scores <- read.xlsx("04_partexdel/splai_part_exdel.xlsx")
exdel_scores <- read.xlsx("04_partexdel/splai_part_exdel_2.xlsx")
exdel_splaiTrack <- DataTrack(
  range = exdel_scores, chromosome = chr, genome = gen, name="SpliceAI ∆Score", 
  background.title = "#F8ACAC", type = "histogram", 
  ylim = c(-1, 1), yTicksAt = c(-1, 0, 1),
  baseline = 0, col.baseline = "#838383", lty.baseline = 1, lwd.baseline = 1.0,
  groups = c("AG", "AL", "DG", "DL", "Var"),
  col = c("#4D4298", "#E60011", "#06A384", "#AB961D", "black"), legend = FALSE, 
  size = 1
)

ov <- OverlayTrack(
  list(dummyTrack, gmtrack, exdel_splaiTrack, exdelTrack), 
  background.title = "#838383"
)
plotTracks(c(axTrack, ov), from = 985, to = 1155)


#5. Partial intron retention
genemodel <- read.xlsx("05_partintret/model_base.xlsx")
gmtrack <- GeneRegionTrack(
  genemodel, fill = "#383838",
  chromosome = chr, genome = "hg19", size = 0.1, alpha = 0.60)

partintret <- read.xlsx("05_partintret/model_partintret.xlsx")
partintretTrack <- GeneRegionTrack(
  partintret, fill = "red",
  chromosome = chr, genome = "hg19", size = 0.3, alpha = 0.2)

partintret_scores <- read.xlsx("05_partintret/splai_part_intret.xlsx")
partintret_splaiTrack <- DataTrack(
  range = partintret_scores, chromosome = chr, genome = gen, name="SpliceAI ∆Score", 
  background.title = "#F8ACAC", type = "histogram", 
  ylim = c(-1, 1), yTicksAt = c(-1, 0, 1),
  baseline = 0, col.baseline = "#838383", lty.baseline = 1, lwd.baseline = 1.0,
  groups = c("AG", "AL", "DG", "DL", "Var"),
  col = c("#4D4298", "#E60011", "#06A384", "#AB961D", "black"), legend = FALSE, 
  size = 1
)

ov <- OverlayTrack(
  list(dummyTrack, gmtrack, partintret_splaiTrack, partintretTrack), 
  background.title = "#838383"
)
plotTracks(c(axTrack, ov), from = 985, to = 1155)

