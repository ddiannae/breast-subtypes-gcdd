library(data.table)
library(tools)

args <- commandArgs(trailingOnly = TRUE)
FILE <- args[1]
NSIZE <- args[2]
SSIZE <- args[3]
#FILE <- "/home/diana/Workspace/rnapipeline/breast-subtypes/data/exp-matrices/LumB_cpm10_arsyn.tsv"
#NSIZE <- 5
#SSIZE <- 105
setwd(dirname(FILE))
exp.matrix <- fread(basename(FILE), header = T)

d <- lapply(1:NSIZE, function(i){
  sel.columns <- c("gene", sample(colnames(exp.matrix)[-1], SSIZE))
  sample.matrix <- exp.matrix[, ..sel.columns]
  fwrite(sample.matrix, file = paste(basename(file_path_sans_ext(FILE)), "-", i, "-",
                                     format(Sys.time(),'%m_%d-%H_%M_%S'), ".txt", sep = ""), row.names = F)
  
})


