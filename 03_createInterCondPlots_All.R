library(ggplot2)
library(reshape2)
library(RColorBrewer)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
start <- Sys.time()
cat("Starting at, ", format(start, "%H:%M:%OS3"), "\n" )

setwd("/labs/csbig/subtipos_mama")

chrs <- c(as.character(1:22), "X")
c1 <- "healthy"
c2 <- "luma"
c3 <- "lumb"
c4 <- "her2"
c5 <- "basal"

g <- parallel::mclapply(X = chrs, mc.cores = 6,  mc.cleanup = FALSE, FUN = function(ch) {
  cat("Working with chromosome ", ch, "\n")
  cat("Loading data\n")
  
  load(file=paste("rdata/", c1, "/inter-norm/", c1, "_norm_ch", ch, ".RData", sep="")) 
  ml1 <- ml
  ml1$condition <- c1
  load(file=paste("rdata/", c2, "/inter-norm/", c2, "_norm_ch", ch, ".RData", sep="")) 
  ml2 <- ml
  ml2$condition <- c2
  load(file=paste("rdata/", c3, "/inter-norm/", c3, "_norm_ch", ch, ".RData", sep="")) 
  ml3 <- ml
  ml3$condition <- c3
  load(file=paste("rdata/", c4, "/inter-norm/", c4, "_norm_ch", ch, ".RData", sep="")) 
  ml4 <- ml
  ml4$condition <- c4
  load(file=paste("rdata/", c5, "/inter-norm/", c5, "_norm_ch", ch, ".RData", sep="")) 
  ml5 <- ml
  ml5$condition <- c5
  
  ml <- rbind(ml1, ml2, ml3, ml4, ml5)
  ml.ml <- melt(ml, value.name = "MI")
  
  png(paste("All_ch", ch, ".png", sep=""), width = 1200, height = 600)
  myplot <- ggplot(data=ml) + 
    geom_density(aes(x = MI, group =interaction(condition, Chr) , colour = condition)) +
    coord_cartesian(xlim = c(0, 0.25)) + scale_color_manual(values= rev(myPalette(5)))
  print(myplot)
  dev.off() 
})

end <- Sys.time()
cat("Whole thing took ", (end - start), ".\n")
cat("Ending at, ", format(end, "%H:%M:%OS3"), "\n" )

    
