library(data.table)
start <- Sys.time()

setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes")

cat("Loading data\n")
load("rdata/annot.RData")

chrs <- c(as.character(1:22), "X")
conds <- c("healthy", "luma", "lumb", "her2", "basal")

getMatricesByCondAndChrs <- function (){
  load("data_subtypes_arsyn.RData")
  for (cond in conds) {
    ## Extract intra-chromosmal interactions
    M <- data.subtypes.arsyn$M[, data.subtypes.arsyn$Targets$Subtype == cond]
    g <- parallel::mclapply(X = chrs, mc.cores = 6, FUN = function(chr) {
      
      cat("Working with chromosome", chr, "\n")
      genes.annot <- annot[annot$Chr == chr, ]
      #genes <- genes.annot[, c("symbol", "start")]
      genes <- genes.annot[, "EnsemblID", drop = FALSE]
      interactions <- merge(genes, M, by = "row.names")
      rownames(interactions) <- interactions$EnsemblID
      interactions <- interactions[, colnames(M)]
      save(interactions, file = paste("intra/", cond, "/", cond,"_chr_", chr, ".RData", sep = ""), compress="xz")
      #interactions <- interactions[, c("symbol", colnames(M), "start")]
      #interactions <- interactions[order(interactions$start), ]
      #write.table(interactions, file = paste("intra/", cond, "/", cond,"_chr_", chr, ".txt", sep=""), 
      #            sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
    })
  }
}

## Execute order: 1
MISifToMIMatrix <- function(){
  genes <- annot$EnsemblID
  ngenes <- length(genes)
  MIunos <- data.table( source = genes, target = genes, MI = rep(1.0, length(genes)))
  for (cond in conds) {
    cat("Working with condition ", cond, "\n")
    MIvals <- fread(paste(cond, "_MI.txt", sep =""), header = F, sep="\t", nThread = 7)
    colnames(MIvals) <- c("source", "target", "MI")
    MIvals <- rbind(MIvals, MIunos)
    MIvals <- MIvals[order(source, target)]
    MIvals <- data.table(matrix(MIvals$MI, nrow = ngenes, ncol = ngenes, byrow = T, dimnames = list(genes, genes)), keep.rownames=T)
    fwrite(MIvals, file = paste("mi-matrices", paste(cond, "MI_matrix.txt", sep = "_"), sep="/"), 
           sep = "\t", row.names = F, col.names = TRUE, quote = F)
    }
}

library(ComplexHeatmap)
library(circlize)
## Execute order: 2
getMIMatricesByChr <- function() {
  library(ComplexHeatmap)
  library(circlize)
  for (cond in conds) {
    cat("Working with condition ", cond, "\n")
    MIvals <- fread( file = paste("mi-matrices", paste(cond, "MI_matrix.txt", sep = "_"), sep="/"), 
                     header = T, sep="\t", nThread = 7, data.table = F)
    rownames(MIvals) <- MIvals$rn
    MIvals <- MIvals[, -c(1)]
    MIvals <- data.matrix(MIvals)
    condvals <- lapply(X = chrs, FUN = function(chr) {
      cat("Working with chromosome", chr, "\n")
      genes.annot <- annot[annot$Chr == chr, ]
      genes.annot <- genes.annot[order(genes.annot$Start), ]
      genes <- genes.annot$EnsemblID
      chr.mi.matrix <- MIvals[genes, genes]
      
      h1 <- Heatmap(chr.mi.matrix, col = colorRamp2(c(0, 1), c("white", "blue")), column_title = paste(cond, "chr", chr),  column_title_gp = gpar(fontsize = 24, fontface = "bold"),
                    cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = F, show_column_names = F, show_row_names = F)
      
      png(paste("mi-matrices/bychr/", cond, "/chr_", chr, "_mi.png", sep=""), width = 1500, height = 1500)
      draw(h1) 
      dev.off()
      write.table(chr.mi.matrix, file =  paste("mi-matrices/bychr/", cond, "/chr_", chr, "_mi_2.txt", sep = ""), 
                  row.names = T, col.names = T, sep = "\t", quote = F)  
  
    })
  }
}

## Execute order: 3
getMIDistanceListByCondAndChrs <- function(){
  cat("getMIDistanceListByCondAndChrs")
  for (cond in conds) {
    cat("Working with condition ", cond, "\n")
    condvals <- lapply(X = chrs, FUN = function(chr) {
      cat("Working with chromosome", chr, "\n")
      MIvals <-read.delim( file = paste("mi-matrices/bychr/", cond, "/chr_", chr, "_mi.txt", sep = ""), 
                           row.names = 1, header = T, sep = "\t")  
      genes <- rownames(MIvals)
      genes.annot <- annot[genes,  ]
      ngenes <- length(genes)
      chrvals <- parallel::mclapply(X = 1:(ngenes-1), mc.cores = 7,  mc.cleanup = FALSE, FUN = function(i) {
        gene1 <- genes[i]
        other.genes <- genes[(i+1):ngenes]
        mivals <- lapply(X = other.genes, 
                         FUN = function(gene2) {
                           distance <- max(genes.annot[gene2, "Start"], genes.annot[gene1, "Start"]) - min(genes.annot[gene1, "Start"], genes.annot[gene2, "Start"])
                           mi <- MIvals[gene1, gene2]
                           data.frame(source = gene1, target = gene2, distance = distance, mi = mi)
                         })
          plyr::ldply(mivals)
      })
      chrdf <- plyr::ldply(chrvals)
      chrdf$chr <- chr
      write.table(chrdf, file =  paste("intra/", cond, "/chr_", chr, "_distance_mi.txt", sep = ""), 
                  row.names = F, col.names = T, sep = "\t", quote = F)  
      chrdf
    })
    conddf <- plyr::ldply(condvals)
    conddf$cond <- cond
    cat("Saving file.\n")
    fwrite(conddf, file =  paste("intra/", cond,"_all_distance_mi.txt", sep = ""), 
                row.names = F, col.names = T, sep = "\t")  
  }
}

getAllMIDistanceMeans <- function(binsize) {
  conditiondfs <- lapply(X = conds, FUN = function(cond){
    cat("Working with condition ", cond, "\n")
    
    meansbych <- parallel::mclapply(X = chrs, mc.cores = 7,  mc.cleanup = FALSE, FUN = function(chr){
      dist.mi.df <- read.delim(file=paste("intra/", cond,"/", "chr_", chr, "_distance_mi.txt", 
                                          sep=""), header = T)
      dist.mi.df <- dist.mi.df[order(dist.mi.df$distance), ]
      rownames(dist.mi.df) <- NULL
      dist.mi.df$bin <- ((as.numeric(rownames(dist.mi.df)) - 1)%/%binsize) + 1
      dfmeans <- aggregate(cbind(distance, mi)~bin, data=dist.mi.df, FUN=mean, na.rm=TRUE)
      dfmeans$cond <- cond
      dfmeans$ch <- chr
      return(dfmeans)
    })
    df <- plyr::ldply(meansbych)
    write.table(df, file = paste("intra/intra-fixed-bin-size/",  binsize, "/", cond, "_all.txt", sep=""), sep="\t", 
                col.names = T, row.names = F, quote = F)
    return(df)
  })
  
  cdf <- plyr::ldply(conditiondfs)
  write.table(cdf, file = paste("intra/intra-fixed-bin-size/", binsize, "/all.txt", sep=""), sep="\t", 
              col.names = T, row.names = F, quote = F)
}

getAllMIDistanceMeans(100)

end <- Sys.time()
cat("Whole thing took ", (end - start), "\n")
cat("Ending at, ", format(end, "%H:%M:%OS3"), "\n" )
