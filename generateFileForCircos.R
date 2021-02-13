library(data.table)

setwd("/mnt/ddisk/transpipeline-data/breast-data/subtypes")

load("rdata/annot.RData")
annot <- annot[!duplicated(annot$EnsemblID), ]
annot <- annot[annot$EnsemblID != "", ]
annot <- annot[, c("EnsemblID", "Chr", "Start", "End")]
conds <- c("LumA", "LumB", "Her2", "Basal", "Healthy")

getMILinks <- function() {
  for(cond in conds) {
    MIvals <- fread(paste("networks/aracne-cluster/1e8/", cond, "-13319.tsv", sep =""),
                    header = T, sep="\t")
    colnames(annot) <-  c("Source", "Source.Chr", "Source.Start", "Source.End")
    MIvals <- merge(annot, MIvals, by.x = "Source", by.y = "Source")
    colnames(annot) <-  c("Target", "Target.Chr", "Target.Start", "Target.End")
    MIvals <- merge(annot, MIvals, by.x = "Target", by.y = "Target")
    MIvals$Source.Chr <- paste0("hs", MIvals$Source.Chr)
    MIvals$Target.Chr <- paste0("hs", MIvals$Target.Chr)
    MIvals <- MIvals[, c("Source.Chr", "Source.Start", "Source.End", 
                         "Target.Chr", "Target.Start", "Target.End")]
    fwrite(MIvals, paste0("circos/", cond, "-circos-data.tsv"),
           sep = "\t", col.names = F,row.names = F)
  }
}

getLFCGraph <- function(){
  for(cond in conds) {
    LFCvals <- fread(paste("deg/treat-", cond, ".tsv", sep =""),
                    header = T, sep="\t", select = c(1, 3, 4), col.names = c("EnsemblID", "lfc", "p"))
    LFCvals <- LFCvals[LFCvals$p <= 0.05 & abs(LFCvals$lfc) > 1, ]
    LFCvals <- merge(annot, LFCvals, by = "EnsemblID")
    LFCvals$Chr <- paste0("hs", LFCvals$Chr)
    LFCvals <- LFCvals[, c("Chr", "Start", "End", "lfc")]
    fwrite(LFCvals, paste0("circos/", cond, "-circos-lfc-filtered-data.tsv"),
           sep = "\t", col.names = F,row.names = F)
  }
}

getExpGraph <- function(){
  for(cond in conds) {
    LFCvals <- fread(paste("exp-matrices/", cond, "_cpm10_arsyn.tsv", sep =""),
                     header = T, sep="\t")
    genes <- LFCvals$gene
    LFCvals <- LFCvals[,-1 ]
    #log2.LFCvals <- log2(LFCvals)
    log2.LFCvals <- LFCvals/10000
    log2.exp.means <- rowMeans(log2.LFCvals)
    log2.exp.means <- data.frame(EnsemblID = genes, Exp = log2.exp.means)
    # log2.max.mean <- max(log2.exp.means$Exp)
    #log2.exp.means$Exp <- log2.exp.means$Exp/log2.max.mean
    LFCvals <- merge(annot, log2.exp.means, by = "EnsemblID")
    LFCvals$Chr <- paste0("hs", LFCvals$Chr)
    LFCvals <- LFCvals[, c("Chr", "Start", "End", "Exp")]
    fwrite(LFCvals, paste0("circos/", tolower(cond), "-circos-exp-data.tsv"),
           sep = "\t", col.names = F,row.names = F)
  }
}

getCNVGraph <- function(){
  for(cond in conds) {
    CNVvals <- fread(paste("cnvs/", tolower(cond), "-gistic/scores.gistic", sep =""),
                     header = T, sep="\t")
    CNVvals <- CNVvals[CNVvals$`-log10(q-value)` > 0, ]
    CNVvals$Chromosome <- paste0("hs", CNVvals$Chromosome)
    CNVvals$`G-score` <- ifelse(CNVvals$Type == "Amp",  CNVvals$`G-score` ,  CNVvals$`G-score` *-1)
    CNVvals <- CNVvals[, c("Chromosome", "Start", "End", "G-score")]
    fwrite(CNVvals, paste0("circos/", tolower(cond), "-circos-cnvs-data.tsv"),
           sep = "\t", col.names = F,row.names = F)
  }
}


