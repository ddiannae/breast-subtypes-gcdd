args = commandArgs(trailingOnly=TRUE)
DIR <- args[1]

files = list.files(path = DIR, pattern = "*.RData", full.names = TRUE, recursive = TRUE)

l <- lapply(files, function(f) {
  cat("Working with file ", f, "\n")
  load(file = f)
  write.table(dfmeans, file = paste(tools::file_path_sans_ext(f),".tsv", sep=""),
                                    quote = F, row.names = F, col.names = T, sep = "\t")
  
})
