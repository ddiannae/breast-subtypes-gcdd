library(dplyr)
chrs <- c(as.character(1:22), "X")
setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes/")

annot <- read.delim("/mnt/ddisk/transpipeline-data/biomarts/Biomart_EnsemblG94_GRCh38_p12.txt", stringsAsFactors = F)
annot <- annot[!duplicated(annot$Gene.stable.ID), ]
annot <- annot[annot$HGNC.symbol != "", ]
colnames(annot) <- c("EnsemblID", "Chr", "Start", "End", "GC", "Type", "Symbol")

genelist <- read.table("networks/genelist.txt", stringsAsFactors = F, col.names = c("EnsemblID"))
genelist <- genelist %>% inner_join(annot, by = "EnsemblID")

N <- nrow(genelist)
N.2 <- choose(N, 2)

all.chrs.p <- lapply(chrs, function(chr){
  xi <- genelist %>% filter(Chr == chr) %>% nrow()    
  return(choose(xi, 2))
})
all.chrs.p <- unlist(all.chrs.p)
names(all.chrs.p) <- chrs
all.chrs.p <- all.chrs.p/N.2


