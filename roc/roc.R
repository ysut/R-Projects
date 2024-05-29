install.packages("pROC")
install.packages("ggplot2")
library(pROC)
library(ggplot2)

input_tsv <- "/Volumes/vol/work/workspace/ROC_Bootstrap_Method/data_3.tsv"
df <- read.table(input_tsv, header = TRUE, sep = "\t")

ROCs <- roc(LABEL ~ PriorityScore + maxsplai, data = df, ci = TRUE, ci.method = "bootstrap")

ROCs

ggroc(ROCs, 
      aes=c("linetype","color"),
      size=1.5,
      legacy.axes = TRUE) +
  geom_abline(color="dark grey", size=0.5)


roc1 <- roc(LABEL~PriorityScore, data = df)
roc2 <- roc(LABEL~maxsplai, data = df)
roc.test(roc1, roc2, method="delong")
roc.test(roc1, roc2, method="bootstrap")


