library(ggplot2)
library(RColorBrewer)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
start <- Sys.time()
cat("Starting at, ", format(start, "%H:%M:%OS3"), "\n" )

setwd("/labs/csbig/subtipos_mama")
#setwd("~/Workspace/subtipos_mama")

chrs <- c(as.character(1:22), "X")
conditions <- c("healthy", "luma", "lumb", "her2", "basal")
l <- lapply(conditions, function(cond) {
  g <- parallel::mclapply(X = chrs, mc.cores =6,  mc.cleanup = FALSE, FUN = function(ch) {
    cat("Working with chromosome ", ch, "\n")
    cat("Loading data\n")
    
    load(file=paste("rdata/", cond, "/inter-take2/", cond, "_ch", ch, ".RData", sep="")) 
    ml
  })
  g <- plyr::ldply(g)
  png(paste("Inter_", cond, ".png", sep=""), width = 1200, height = 600)
  myplot <- ggplot(data=g) + 
      geom_density(aes(x = MI, colour = myPalette(1)), show.legend = F) +
      coord_cartesian(xlim = c(0, 0.25)) 
  print(myplot)
  dev.off() 
})
end <- Sys.time()
cat("Whole thing took ", (end - start), ".\n")
cat("Ending at, ", format(end, "%H:%M:%OS3"), "\n" )

    
