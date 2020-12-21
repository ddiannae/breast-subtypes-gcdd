start <- Sys.time()
cat("Starting at, ", format(start, "%H:%M:%OS3"), "\n" )

setwd("~/Workspace/subtipos_mama")
#setwd("/labs/csbig/subtipos_mama")

cat("Loading data\n")
load(file="genes.annot.RData")
load(file="adjmatrix_basal.RData")

all.genes.annot <- genes.annot
chrs <- c(as.character(1:22), "X")
binlength <- 100000

## Extract intra-chromosmal interactions
for (ch in chrs) {
  cat("Working with chromosome", ch, "\n")
  genes.annot <- all.genes.annot[all.genes.annot$Chr == ch, ]
  genes <- genes.annot$symbol
  g <- parallel::mclapply(X = genes, mc.cores = 6,  mc.cleanup = FALSE, FUN = function(gene1) {
    #cat("Gene: ", gene1, "\n")
    other.genes <- genes[which(genes != gene1)]
    
    mivals <- lapply(X = other.genes, 
                                FUN = function(gene2) {
                                  distance <- max(genes.annot[gene2, "start"], genes.annot[gene1, "start"]) - min(genes.annot[gene1, "end"], genes.annot[gene2, "end"])
                                  bin <- max(floor(distance/binlength), 1)
                                  list(bin, max(adjmatrix[gene1, gene2], adjmatrix[gene2, gene1], na.rm = T))
                                })
    matrix(unlist(mivals), ncol = 2, byrow = T)
  })
  cat("Getting bin order.\n")
  g <- plyr::ldply(g)
  colnames(g) <- c("bin", "mi")
  # 
  # o <- lapply(unique(g$bin), function(b) {
  #   l <- g[g$bin == b, ]
  #   l$order <- order(l$mi)
  #   l
  # })
  # 
  # o <- plyr::ldply(o)
  # o <- o[order(o$bin, o$order),  ]
  cat("Saving file.\n")
  write.table(g, file=paste("lumb_intra_order_",ch , ".tsv", sep=""), row.names = F, col.names = T, sep = "\t")
}
end <- Sys.time()
cat("Whole thing took ", (end - start), "\n")
cat("Ending at, ", format(end, "%H:%M:%OS3"), "\n" )

