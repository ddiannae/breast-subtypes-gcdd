library(ggplot2)
start <- Sys.time()
cat("Starting at, ", format(start, "%H:%M:%OS3"), "\n" )

setwd("~/Workspace/subtipos_mama/")

cat("Loading data\n")
load(file="rdata/lista_distancias_intracromosomales.RData")
conds <- c("basal", "healthy")
f <- 0.001
ff <- "001"

for (cond in conds) {
  
  cat("Working with condition ", cond, "\n")
  load(file=paste("rdata/adjmatrix_", cond, ".RData", sep=""))
  
  chrs <- c(as.character(1:22), "X")
  binsizes <- c() 
  for (ch in chrs) {
    cat("\t Working with chromosome ", ch, "\n")
    distancias <- lista_distmatrices[[ch]]
    genes <- rownames(distancias)
    r <- rep(genes, each = length(genes))
    c <- rep(genes, times = length(genes))
    dist.df <- data.frame(source = r, target = c, distance = unlist(distancias),  row.names = NULL)
    dist.df <- dist.df[!is.na(dist.df$distance), ]
    dist.df <- dist.df[dist.df$distance > 0, ]
    
    mis <- adjmatrix[rownames(adjmatrix) %in% genes, colnames(adjmatrix) %in% genes]
    genes <- rownames(mis)
    r <- rep(genes, each = length(genes))
    c <- rep(genes, times = length(genes))
    mi.df <- data.frame(source = r, target = c, mi = unlist(mis), row.names = NULL)
    mi.df <- mi.df[!is.na(mi.df$mi), ]
    
    dist.mi.df <- merge(dist.df, mi.df)
    dist.mi.df <- dist.mi.df[with(dist.mi.df, order(distance, mi)), ]
    
    cat("Valor máximo duplicados: ", max(table(dist.mi.df[duplicated(dist.mi.df$distance), "distance"])), "\n")
    
    #write.table(dist.mi.df, file = paste(cond,"_mi_distance_", ch, ".tsv", sep=""), sep="\t", 
    #            col.names = T, row.names = F, quote = F)
    
    binsize <- floor(dim(dist.mi.df)[1]*f)
    binsizes <- c(binsizes, binsize)
    bins <- split(dist.mi.df, seq(from =1, to = nrow(dist.mi.df)-1)%/%binsize)
    testy <- parallel::mclapply(X = seq(1:length(bins)), mc.cores = 6,  mc.cleanup = FALSE, FUN = function(i){
      n <- names(bins[i])
      b <- bins[[i]]
      #png(paste("plots/intra-fixed-bins/binsdensity/", cond, "_ch", ch, "_bin", n, ".png", sep=""), width =800 , height = 400)
      myplot <- ggplot(b, aes(mi))+geom_density()
      p <- layer_data(myplot)
      # myplot <- myplot +
      #   geom_vline(aes(xintercept=mean(b$mi)), color="blue", size=1) +
      #   geom_vline(aes(xintercept=mean(p[which.max(p$y), "x"])), color="red", size=1)
      # print(myplot)
      # dev.off()
      return(matrix(c(p[which.max(p$y), "x"], mean(b$mi), median(b$mi), max(b$mi)), ncol = 4))
    })

    stats <- plyr::ldply(testy)
    colnames(stats) <- c("moda", "media", "mediana", "max")
    write.table(stats, file = paste("rdata/intra-fixed-bins/stats_", cond,"_ch", ch, "_", binsize, ".tsv", sep=""), sep="\t",
                col.names = T, row.names = T, quote = F)
    #stats <- read.delim(file=paste("rdata/intra-fixed-bins/stats_", cond,"_", binsize, "_ch", ch, ".tsv", sep=""),
    #                    header = T)
    s <- "media"
    #for (s in colnames(stats)) {
    png(paste("plots/intra-fixed-bins/", cond, "_ch", ch, "_", binsize, "_", s, ".png", sep=""), width =800 , height = 400)
    myplot <- ggplot(data=stats, aes(x=as.numeric(rownames(stats)), y=stats[,s], group=1)) +
      geom_point(size = 1.5) + labs(x = "B in", y = s)
    print(myplot)
    dev.off()
    #}
    
  }
  # write.table(binsizes, file = paste("rdata/intra-fixed-bins/binsizes_", ff,".tsv", sep=""), sep="\t",
  #         row.names = T, quote = F)
}

