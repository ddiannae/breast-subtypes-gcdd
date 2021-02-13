library(data.table)
setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes/")

conds <- c("healthy", "luma", "lumb", "basal", "her2")

annot <- read.delim("/mnt/ddisk/transpipeline-data/biomarts/Biomart_EnsemblG94_GRCh38_p12.txt", stringsAsFactors = F, 
                    col.names = c("EnsemblID", "Chr", "Start", "End", "GC",  "Type", "Symbol") )
annot <- annot[!duplicated(annot$EnsemblID), ]
chrs <- c(as.character(1:22), "X")
conds <- c("healthy", "luma", "lumb", "her2", "basal")

for (cond in conds) {
  cat("Working with condition ", cond, "\n")
  MImatrix <-fread(file = paste0("networks/aracne-cluster/1/", cond, ".tsv", sep = ""))
  MImatrix <- data.matrix(MImatrix)
  diag(MImatrix) <- 0
  maxMI <- max(MImatrix, na.rm = T)
  MImatrix <- MImatrix/maxMI
  diag(MImatrix) <- 1
  fwrite(MImatrix, file =  paste0("networks/aracne-cluster/1/", cond, "_normMI.tsv"),
         row.names = F, col.names = T, sep = "\t")  
}
