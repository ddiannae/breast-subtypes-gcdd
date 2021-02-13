library(data.table)
suppressMessages(library("car"))
suppressMessages(library("ggplot2"))
suppressMessages(library("reshape"))
suppressMessages(library("cowplot"))
suppressMessages(library("grid"))

setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes/")


load("rdata/annot.RData")
annot <- annot[!duplicated(annot$EnsemblID), ]
annot <- annot[annot$EnsemblID != "", ]
annot <- annot[, c("EnsemblID", "Chr", "Start", "End")]
conds <- c("LumA", "LumB", "Her2", "Basal")


LFCall <- lapply(conds, function(cond){
  LFCvals <- fread(paste("deg/ebayes-", cond, ".tsv", sep =""),
                   header = T, sep="\t", select = c(1, 3, 4), col.names = c("EnsemblID", "lfc", "p"))
  LFCvals <- LFCvals[LFCvals$p <= 0.05 & abs(LFCvals$lfc) > 1, ]
  LFCvals <- merge(annot, LFCvals, by = "EnsemblID")
  LFCvals <- LFCvals[, c("Chr", "Start", "lfc")]
  LFCvals$Subtype <- cond
  return(LFCvals)
})
LFCall <- rbindlist(LFCall)
LFCall$col <- ifelse(LFCall$lfc < 0, "down", "up")
LFCall$col <- as.factor(LFCall$col)
basal <- LFCall[LFCall$Subtype == "Basal", ]
basal.1 <- basal[basal$Chr == "1", ]
basal.8 <- basal[basal$Chr == "8", ]
basal.8 <- basal.8[order(-basal.8$lfc), ]
basal.1.pos <- basal.1[basal.1$lfc > 0, ]

LFCall$Chr <- factor(LFCall$Chr, levels = c(as.character(1:22), "X"))
LFCall$Subtype <-  factor(LFCall$Subtype, levels = c("Healthy", "LumA", "LumB", "Her2", "Basal"))

p <- ggplot(
  data = LFCall,
  aes(
    x = Start,
    y = lfc,
    color = col
  )
) + 
  geom_bar(stat = "identity", position = "identity", lwd = 0.5) +
  scale_color_manual(values=c(down="blue", up="red")) + 
  facet_grid(Chr ~ Subtype) +
  theme_bw(base_size = 20) +
  theme(
      legend.position = "none",
      strip.text.x = element_text(size = 30),
      strip.text.y  = element_text(size = 24),
      axis.title=element_text(size=30,face="bold")
    ) + xlab("Start") + ylab("LFC")





+


+
  guides(
    color = guide_legend( 
      nrow = 2,
      byrow = TRUE
    )
  )



+
  