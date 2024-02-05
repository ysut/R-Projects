install.packages("GGally")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("ggforce")

# 必要なパッケージを読み込む
library(ggplot2)
library(dplyr)
library(GGally)
library(ggforce)

# サンプルデータの生成
set.seed(123) 
samples <- data.frame(
  ID = 1:1000,
  Status = rep(c("Resolved", "Unresolved"), c(400, 600)),
  Pathogenicity = c(rep(c("Yes", "No"), c(350, 50)), rep(c("Yes", "No"), c(200, 400)))
)


ggplot(samples, aes(color = Pathogenicity)) +
  geom_parallel_sets(aes(x = Status, id = ID), alpha = 0.5) +
  scale_fill_manual(values = c("Yes" = "blue", "No" = "red")) +
  theme_minimal() +
  labs(title = "Parallel Plot of Sample Status and Pathogenicity",
       x = "Status",
       y = "Frequency",
       color = "Pathogenicity")

