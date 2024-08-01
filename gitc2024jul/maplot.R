
arab2 <- "/Volumes/vol/gitc/data/SS/arab2.txt"
dat <- read.delim(arab2, row.names=1)

colSums(dat)

plot(dat$m1 + 1, dat$m2 + 1, log="xy")
dat$m1

M <- log2(dat$m2 +1 ) - log2(dat$m1 +1)
A <- 1/2 * (log2(dat$m2 + 1) + log2(dat$m1 + 1))
plot(A, M)
plot(A, M, pch=16, cex=0.4, ylim=c(-8,8), main="MA plot (m1 vs m2)")
abline(h=log2(2), col="red", lty=2)
abline(h=log2(1/2), col="red", lty=2)
abline(h=0, col="blue", lty=2)
