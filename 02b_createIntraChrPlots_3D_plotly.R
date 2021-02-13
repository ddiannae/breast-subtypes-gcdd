library(plotly)
library(tibble)
library(dplyr)
Sys.setenv("plotly_username"="dianagco")
Sys.setenv("plotly_api_key"="jh47zwssSzbvn6gUsfAE")

start <- Sys.time()
cat("Starting at, ", format(start, "%H:%M:%OS3"), "\n" )

setwd("~/Workspace/subtipos_mama")

#conditions <- c("healthy", "luma", "lumb", "her2", "basal")
cond <- "healthy"
ch <- as.character(1)

mi.df <- read.delim(paste(cond, "_intra_order_", ch, ".tsv", sep=""),
                        stringsAsFactors = T, header = T, sep="\t" )
  
mi.df <- mi.df[order(mi.df$bin, -mi.df$mi),]
mi.df$orden <- NULL
  
aux <- tapply(
  mi.df$bin, 
  INDEX = factor(mi.df$bin, levels = 1:max(mi.df$bin)),
  function(x){
    1:length(x)
  }
)
aux <- do.call(c, aux)
mi.df$orden <- aux
mi.df <- mi.df[c(1,3,2)]
mi.agg <- aggregate(mi ~ orden, data=mi.df, FUN = cbind)
## max_length <- max(unlist(lapply (mi.agg$mi, FUN = length)))
## Solo tomemos los primeros 2000 por bin 
max_length <- 400
mi.matrix <- sapply (mi.agg$mi, function (x) {
                length (x) <- max_length; 
                return (x)
              })
#mi.matrix[is.na(mi.matrix)] <- 0 
#min.matrix <- mi.matrix[1:1500, 1:1000]
min.matrix <- mi.matrix

ax <- list(
  zeroline = FALSE,
  showline = T,
  showticklabels = FALSE,
  showgrid = T
)

p <- plot_ly(z = ~min.matrix, type = 'surface',
              colorbar=list(title=''),
              colorscale =  list(c(0, 0.05, 0.1, 0.2, 0.5, 1), 
                                 c("rgb(4, 35, 51)", "rgb(60, 56, 143)",
                                   "rgb(137, 82, 141)", "rgb(212, 108, 107)",
                                   "rgb(249, 165, 66)", "rgb(232, 250, 91)"))) %>%
    layout(scene = list(xaxis=c(ax, title="Order"), 
                                yaxis = c(ax, title = "Distance"), 
                                zaxis= c(ax, title = "MI value", showline = F)))

p

