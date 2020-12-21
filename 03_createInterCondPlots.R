library(ggplot2)
library(reshape2)
library(RColorBrewer)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
start <- Sys.time()
cat("Starting at, ", format(start, "%H:%M:%OS3"), "\n" )

setwd("/labs/csbig/subtipos_mama")

chrs <- c(as.character(1:22), "X")
c1 <- "healthy"
c2 <- "normal"

g <- parallel::mclapply(X = chrs, mc.cores = 5,  mc.cleanup = FALSE, FUN = function(ch) {
  cat("Working with chromosome ", ch, "\n")
  cat("Loading data\n")
  
  load(file=paste("rdata/", c1, "/inter-take2/", c1, "_ch", ch, ".RData", sep="")) 
  ml1 <- ml
  ml1$condition <- c1
  load(file=paste("rdata/", c2, "/inter-take2/", c2, "_ch", ch, ".RData", sep="")) 
  ml2 <- ml
  ml2$condition <- c2
  
  ml <- rbind(ml1, ml2)
  ml.ml <- melt(ml, value.name = "MI")
  
  png(paste(c1, "_", c2, "_ch", ch, ".png", sep=""), width = 1200, height = 600)
  myplot <- ggplot(data=ml) + 
    geom_density(aes(x = MI, group =interaction(condition, Chr) , colour = Chr)) +
    coord_cartesian(xlim = c(0, 0.25)) + scale_color_manual(values= rev(myPalette(23)))
  print(myplot)
  dev.off() 
})

end <- Sys.time()
cat("Whole thing took ", (end - start), ".\n")
cat("Ending at, ", format(end, "%H:%M:%OS3"), "\n" )

    
