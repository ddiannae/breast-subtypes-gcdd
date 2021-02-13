setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes/networks/aracne-cluster/1e8/13319")

conds <- c("basal", "luma", "lumb", "her2", "healthy")

for (cond in conds) {
  network <- read.delim(paste0(cond,"-13319.tsv"), stringsAsFactors = F)
  maxMI <- max(network$MI)
  network$MI <- network$MI / maxMI
  write.table(network, paste0(cond,"-13319-norm-mi.tsv"), sep = "\t", quote = F, row.names = F)
}

