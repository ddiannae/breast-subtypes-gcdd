library(circlize)
library(ComplexHeatmap)

annot <- read.delim("/mnt/ddisk/transpipeline-data/biomarts/Biomart_EnsemblG94_GRCh38_p12_NCBI.txt", stringsAsFactors = F)
annot <- annot[!duplicated(annot$Gene.stable.ID), ]
annot <- annot[annot$HGNC.symbol != "", ]
colnames(annot) <- c("EnsemblID", "Chr", "Start", "End", "GC", "Type", "Symbol", "NCBI")
chrs <- c(as.character(1))
conds <- c("Healthy", "LumA", "LumB", "Her2", "Basal")

## Execute order: 2
getMIMatricesByChr <- function() {
  library(ComplexHeatmap)
  library(circlize)
  for (cond in conds) {
    cat("Working with condition ", cond, "\n")
    MImatrix <-fread( file = paste0("networks/aracne-cluster/1/", cond, ".tsv", sep = ""),
                      header = T, sep = "\t")
    MImatrix <- data.matrix(MImatrix)
    MImatrix[is.na(MImatrix)] <- 0
    MImatrixT <- t(MImatrix)
    diag(MImatrixT) <- 0
    MImatrixD <- MImatrix + MImatrixT

    rownames(MImatrix) = colnames(MImatrix)
    condvals <- lapply(X = chrs, FUN = function(chr) {
      cat("Working with chromosome", chr, "\n")
      genes.annot <- annot[annot$Chr == chr, ]
      genes.annot <- genes.annot[order(genes.annot$Start), ]
      genes <- genes.annot$EnsemblID
      chr.mi.matrix <- MImatrix[genes, genes]

      h1 <- Heatmap(chr.mi.matrix, col = colorRamp2(c(0, 0.2), c("white", "blue")), column_title = paste(cond, "chr", chr),  column_title_gp = gpar(fontsize = 24, fontface = "bold"),
                    cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = F, show_column_names = F, show_row_names = F)

      png(paste("mi-matrices/bychr/", cond, "/chr_", chr, "_mi.png", sep=""), width = 1500, height = 1500)
      draw(h1)
      dev.off()
      write.table(chr.mi.matrix, file =  paste("mi-matrices/bychr/", cond, "/chr_", chr, "_mi_2.txt", sep = ""),
                  row.names = T, col.names = T, sep = "\t", quote = F)

    })
  }
}
