library(ggplot2)
library(cowplot)
setwd("~/Workspace/subtipos_mama/")

conds <-  c("basal", "luma", "lumb", "her2", "healthy")
#conds <-  c("basal")
binsize <- 100

t <- parallel::mclapply(X = conds, mc.cores = 1,  mc.cleanup = FALSE, FUN = function(cond){
  datos_chr1 <- read.delim(file=paste("rdata/pares-distancia-mi/", cond,"/", cond, "_distance_mi1.tsv", sep=""), header = T)
  datos_all <- read.delim(file=paste("rdata/pares-distancia-mi/", cond,"/", cond,"_all.tsv", sep=""), header = T)
  
  #this should be ok
  datos_chr1 <- datos_chr1[with(datos_chr1, order(distance)), ]
  #this is neccesary
  datos_all <- datos_all[with(datos_all, order(distance)), ]
  rownames(datos_all) <- 1:nrow(datos_all)
  
  datos_chr1$bin <- ((as.numeric(rownames(datos_chr1)) - 1)%/%binsize) + 1
  datos_all$bin <- ((as.numeric(rownames(datos_all)) - 1)%/%binsize) + 1
  
  dfmeans_chr1 <- aggregate(mi ~ bin, data=datos_chr1, FUN=mean, na.rm=TRUE)
  dfmeans_all <- aggregate(mi ~ bin, data=datos_all, FUN=mean, na.rm=TRUE)
  
  p_all <- ggplot() +
    geom_point(
      data = dfmeans_all,
      aes(
        x = bin,
        y = mi
      ),
      color = "gray",
      size = 2
    ) +
    geom_smooth(
      data = dfmeans_all,
      aes(
        x = bin,
        y = 10^mi
      ),
      method = 'lm',
      formula = log10(y) ~ 1 + log10(x) + I(log10(x)^2)  + I(log10(x)^3),
      color = "black",
      size=1
    ) +
    xlab("Distance (a. u.)") +
    ylab("Mutual Information")  + 
    theme_classic(base_size = 10)
  
  p_chr1 <- ggplot() +
    geom_point(
      data = dfmeans_chr1,
      aes(
        x = bin,
        y = mi
      ),
      color = "gray",
      size = 2
    ) +
    geom_smooth(
      data = dfmeans_chr1,
      aes(
        x = bin,
        y = 10^mi
      ),
      method = 'lm',
      formula = log10(y) ~ 1 + log10(x) + I(log10(x)^2)  + I(log10(x)^3),
      color = "black", 
      size= 1
    ) +
    xlab("Distance (a. u.)") +
    ylab("Mutual Information")  + 
    theme_classic(base_size = 10)
  #p_chr1
  p <- plot_grid(
    plotlist = list(p_all, p_chr1),
    nrow = 2, ncol = 1
  )
  png(paste("plots/intra-fixed-bins-size/", cond, "_supplementary.png", sep=""), width = 400 , height = 400)
  print(p)
  dev.off()
})
