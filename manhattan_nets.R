library(tidyverse)
source("/home/diana/Workspace/breast-cancer-subtypes/manhattan.R")

setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes")


load("rdata/annot.RData")
chrs <- c(as.character(1:22), "X")
conds <- c("healthy", "luma", "lumb", "her2", "basal")

for (cond in conds) {
  network <- read_tsv(paste0("networks/aracne-cluster/1e8/13319/", cond, "-13319.tsv"))
  network <- select(network, Source, MI, Target)
  network2 <- select(network, Target, MI, Source)
  colnames(network2) <- c("Source", "MI", "Target")
  network <- bind_rows(network, network2)
  rm(network2)
  mi.by.source <- network %>% group_by(Source) %>% summarise(MI.total = sum(MI), Avg.MI = mean(MI), grado = n())
  
  exp.diff <- read_tsv(paste0("deg/treat-", cond, ".tsv"))
  exp.diff <- exp.diff[, c(1,3)]
  colnames(exp.diff) <- c("EnsemblID", "LFC")
  
  all.mi.by.source <- lapply(X = chrs, FUN = function(chr) {
    cat("Working with chromosome", chr, "\n")
    
    genes.annot <- as_tibble(annot[annot$Chr == chr, c("EnsemblID", "Start") ])
    chr.mi.by.source <- inner_join(genes.annot, mi.by.source, by = c("EnsemblID" = "Source"))
    chr.mi.by.source <- inner_join(chr.mi.by.source, exp.diff, by = c("EnsemblID" = "EnsemblID"))
    chr.mi.by.source <- chr.mi.by.source %>% mutate(chr = chr) %>% arrange(Start)
    return(chr.mi.by.source)
  })
  
  all.mi.by.source <- plyr::ldply(all.mi.by.source)
  all.chrs <- all.mi.by.source$chr
  all.chrs <- ordered(all.chrs, levels = chrs)
  
  chr1 <- all.mi.by.source[all.mi.by.source$chr == "1", ]
  
  png(paste0(cond, "_manhatans.png"), width = 2000, height = 900)
  par(mfcol=c(1,3))
  
  manhattan.plot(all.chrs, all.mi.by.source$Start, all.mi.by.source$LFC, ylab = "LFC", col=c("darkorange", "darkgray"), cex = 0.4)
  manhattan.plot(chr1$chr, chr1$Start, chr1$LFC, ylab = "LFC", col=c("darkorange", "darkgray"), cex = 0.4)
  manhattan.plot(all.chrs, all.mi.by.source$Start, all.mi.by.source$grado, ylab = "Grado", col=c("darkorange", "darkgray"), cex = 0.4)
  manhattan.plot(all.chrs, all.mi.by.source$Start, all.mi.by.source$Avg.MI, ylab = "Avg.MI", col=c("darkorange", "darkgray"), cex = 0.4)

  dev.off()
  
  
}
