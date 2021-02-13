suppressMessages(library("car"))
suppressMessages(library("ggplot2"))
suppressMessages(library("reshape"))
suppressMessages(library("cowplot"))
suppressMessages(library("grid")) 
library(data.table)
library(dplyr)

setwd("/mnt/antares/mnt/ddisk/transpipeline-data/breast-data/subtypes/")

load("rdata/annot.RData")
annot <- annot[!duplicated(annot$EnsemblID), ]
annot <- annot[annot$EnsemblID != "", ]
annot <- annot[, c("EnsemblID", "Chr", "Start", "End")]
conds <- c("luma", "lumb", "her2", "basal")


gistic.vals.all <- lapply(conds, function(cond){
  gistic.del <- fread(paste("cnvs/", cond, "-gistic/del_genes.conf_99.txt", sep =""),
                   header = F, sep="\t", nrows = 4)
  gistic.del <- transpose(gistic.del[, -1])
  colnames(gistic.del) <- c("cytoband", "q", "resq", "widepeak")
  gistic.del$q <- as.numeric( gistic.del$q)
  gistic.del$resq <- as.numeric( gistic.del$resq)
  gistic.del$Subtype <- cond
  gistic.del$Type <- "Loss"
  
  gistic.amp <- fread(paste("cnvs/", cond, "-gistic/amp_genes.conf_99.txt", sep =""),
                      header = F, sep="\t", nrows = 4)
  gistic.amp <- transpose(gistic.amp[, -1])
  colnames(gistic.amp) <- c("cytoband", "q", "resq", "widepeak")
  gistic.amp$q <- as.numeric( gistic.amp$q)
  gistic.amp$resq <- as.numeric( gistic.amp$resq)
  gistic.amp$Subtype <- cond
  gistic.amp$Type <- "Gain"
  
  return(rbind(gistic.del, gistic.amp))
})

gistic.vals.all <- rbindlist(gistic.vals.all)
gistic.vals.all <- gistic.vals.all[!is.na(gistic.vals.all$cytoband), ]
#gistic.vals.all$Chromosome <- factor(gistic.vals.all$Chromosome, levels = c(as.character(1:22), "X"))
#gistic.vals.all$Subtype <-  factor(gistic.vals.all$Subtype, levels = c("healthy", "luma", "lumb", "her2", "basal"))

p <- ggplot(
  data = gistic.vals.all,
  aes(
    x = Start,
    y = `G-score`,
    color = Type
  )
) + 
  geom_bar(stat = "identity", position = "identity", lwd = 0.5) +
  scale_color_manual(values=c(Del="blue", Amp="red")) + 
  facet_grid(Chromosome ~ Subtype) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 30),
    strip.text.y  = element_text(size = 24),
    axis.title=element_text(size=30,face="bold")
  ) + xlab("Start") + ylab("LFC")

paper.results <- read.csv(file = "cnvs/13058_2009_2532_MOESM7_ESM.csv", sep = "\t", header = T, stringsAsFactors = F)
basal.paper.results <- paper.results[paper.results$p.value_.basal < 0.05, ]
basal.own <- gistic.vals.all[gistic.vals.all$Subtype == "basal", ]
her2.own <- gistic.vals.all[gistic.vals.all$Subtype == "her2", ]
merge(her2.own, basal.paper.results, by.x = "cytoband", by.y = "Cytoband")
