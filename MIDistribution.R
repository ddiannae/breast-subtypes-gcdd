library(data.table)
library(ggplot2)
setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes/networks/aracne-cluster")

conds <- c("healthy", "her2", "luma", "lumb", "basal")
networkFiles <- paste(conds, "txt", sep = ".")
mitables <- lapply(X = conds, FUN = function(network) {
  cat("Reading network ", network, "\n")
  MIvals <- fread(paste(network, "sif", sep="."), header = F, sep="\t", nThread = 1, col.names =  c("Source", "MI", "Target"))
  MIvals$Cond <- network
  MIvals
})

DT <- rbindlist(mitables)

p <- ggplot(DT, aes(x = MI, color = Cond)) + geom_density()
png("mi-distribution.png", width = 1500, height = 1500)
print(p)
dev.off()

