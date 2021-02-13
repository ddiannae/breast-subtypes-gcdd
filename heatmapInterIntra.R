library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(ComplexHeatmap)
library(circlize)
setwd("/media/ddisk/transpipeline-data/breast-data/subtypes")

chrs <- c(as.character(1:22), "X")
conds <- c("healthy", "luma", "lumb", "her2", "basal")

all.networks <- lapply(conds, function(cond) {
  
  interactions <- read.delim(paste0("networks/network-tables/", cond,"-norm-mi-interactions.tsv"), 
                        stringsAsFactors = F)
  interactionsT <- interactions
  colnames(interactionsT) <- c("Target", "Source", "MI", "InteractionType", "RowNumber")
  interactions <- rbind(interactions, interactionsT)
  interactions <- interactions %>% arrange(RowNumber)
  
  vertices <- read.delim(paste0("networks/network-tables/", cond,"-vertices.tsv"), 
                         stringsAsFactors = F)
  interactions <- interactions %>% 
    inner_join(vertices, by = c("Source" = "EnsemblID")) %>%  
    inner_join(vertices, by = c("Target" = "EnsemblID"), suffix = c(".Source", ".Target")) 
  
  inter.by.chr <- interactions %>% group_by(Chr.Source, Chr.Target) %>% tally()
  inter.matrix <- inter.by.chr %>% spread(Chr.Target , n,  fill = 0) %>% 
    column_to_rownames(var = "Chr.Source")  %>% select(chrs) 
  max.inter <- inter.by.chr %>% filter(Chr.Source != Chr.Target)
  max.intra <- inter.by.chr %>% filter(Chr.Source == Chr.Target)
  return(list(dm = data.matrix(inter.matrix[chrs, ]), 
              max.inter = max(max.inter$n), 
              max.intra = max(max.intra$n), 
              cond = cond))
})
names(all.networks) <- unlist(lapply(all.networks, "[[", "cond"))
max.inter <- max(unlist(lapply(all.networks, "[[", "max.inter")))
max.intra <- max(unlist(lapply(all.networks, "[[", "max.intra")))

intra_col_fun = colorRamp2(c(0, max(max.intra)), c("white", "#2d20a5"))
inter_col_fun = colorRamp2(c(0, max(max.inter)), c("white", "#ff8833"))

lapply(names(all.networks), function(cond) {
  png(paste0("figures/heatmap_", cond, ".png"), width = 2000, height = 2000, units = "px")
  draw(Heatmap(all.networks[[cond]]$dm, inter_col_fun,
          cluster_rows = FALSE, cluster_columns = F, 
          show_heatmap_legend = F, show_column_names = T,
          column_names_side = "top", column_names_rot = -45,
          row_names_rot = -45,
          column_names_centered = TRUE,
          column_labels = paste0(colnames(all.networks[[cond]]$dm), " "),
          row_labels = paste0("   ", rownames(all.networks[[cond]]$dm)),
          row_names_gp = gpar(fontsize = 40), column_names_gp = gpar(fontsize = 40),
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(i == j) {
              grid.rect(x = x, y = y, width = width, height = height, 
                        gp = gpar(col = NA, fill = intra_col_fun(all.networks[[cond]]$dm[i, j])))
              grid.text(sprintf("%d", all.networks[[cond]]$dm[i, j]), x, y, gp = gpar(fontsize = 25), rot = -45)
            } else if(i > j) {
              grid.rect(x = x, y = y, width = width, height = height, 
                        gp = gpar(col = "white", fill = "white"))
            } else {
              grid.text(sprintf("%d", all.networks[[cond]]$dm[i, j]), x, y, gp = gpar(fontsize = 25), rot = -45)
            } 
          }
  ))
  dev.off()
})

lgd = Legend(col_fun = inter_col_fun, title = "-trans",  legend_height = unit(10, "cm"), 
             title_position = "lefttop-rot", labels_gp = gpar(col = "black", fontsize = 20),
             title_gp = gpar(fontsize = 35))
lgd2 = Legend(col_fun = intra_col_fun, title = "-cis",  legend_height = unit(10, "cm"), 
             title_position = "lefttop-rot", labels_gp = gpar(col = "black", fontsize = 20),
             title_gp = gpar(fontsize = 35))
pd = packLegend(lgd, lgd2,  direction = "horizontal", column_gap = unit(5, "mm"))
dev.off()
draw(pd)

# ### Old Heatmpap inter/intra
# all.networks <- rbindlist(all.networks)
# all.inter <- all.networks %>% filter(Inter == T)
# all.intra <- all.networks %>% filter(Inter == F)
# 
# intra.by.chr <- all.intra %>% group_by(Cond, Source.Chr) %>% tally()
# intra.matrix <- intra.by.chr %>% spread(Cond, n)
# intra.matrix$Source.Chr <- factor(intra.matrix$Source.Chr, levels = chrs)
# intra.matrix <- intra.matrix %>% arrange(Source.Chr) %>% select(Source.Chr, conds)
# rnames <- unlist(intra.matrix[, 1])
# intra.matrix <- as.matrix(intra.matrix[, -1])
# rownames(intra.matrix) <- rnames
# cn <- c("Healthy", "LumA", "LumB", "Her2", "Basal")
# Heatmap(intra.matrix,  colorRamp2(c(0, 300), c("white", "#2d20a5")),
#               cluster_rows = FALSE, cluster_columns = F, show_heatmap_legend = T, show_column_names = F,
#               row_names_gp = gpar(fontsize = 25),
#         heatmap_legend_param = list(title = "Cis", title_gp = gpar(fontsize = 20), 
#                                     labels_gp = gpar(fontsize = 18),  legend_height = unit(5, "cm"),
#                                     grid_height = unit(10, "mm")),
#         
#         bottom_annotation = HeatmapAnnotation(
#             Fenotype = cn,
#             col = list(Fenotype = c("Healthy" = "#C7C7C7", "LumA" = "#FF9E4A", 
#                                "LumB" = "#67BF5C", "Her2" = "#729ECE", "Basal" = "#ED665D")),
#             annotation_height = 70, show_legend = F, annotation_name_gp = gpar(fontsize = 20)
#         )
#             )
# 
# inter.by.s.chr <- all.inter %>% group_by(Cond, Source.Chr) %>% tally()
# colnames(inter.by.s.chr) <- c("Cond", "Chr", "n")
# inter.by.t.chr <- all.inter %>% group_by(Cond, Target.Chr) %>% tally()
# colnames(inter.by.t.chr) <- c("Cond", "Chr", "n")
# inter.matrix <- rbind(inter.by.s.chr, inter.by.t.chr)
# inter.matrix <- inter.matrix %>% group_by(Chr, Cond) %>% tally()
# inter.matrix <- inter.matrix %>% spread(Cond, n)
# inter.matrix$Chr <- factor(inter.matrix$Chr, levels = chrs)
# inter.matrix <- inter.matrix %>% replace(is.na(.), 0) %>% arrange(Chr) %>% select(Chr, conds)
# rnames <- unlist(inter.matrix[, 1])
# inter.matrix <- as.matrix(inter.matrix[, -1])
# rownames(inter.matrix) <- rnames
# Heatmap(inter.matrix,  colorRamp2(c(0, 200), c("white", "#ff8833")),  column_title_gp = gpar(fontsize = 24, fontface = "bold"),
#         cluster_rows = FALSE, cluster_columns = F, show_heatmap_legend = T, show_column_names = F,
#         heatmap_legend_param = list(title = "Trans", title_gp = gpar(fontsize = 20), 
#                                     labels_gp = gpar(fontsize = 18),  legend_height = unit(5, "cm"),
#                                     grid_height = unit(10, "mm")),
#         
#         bottom_annotation = HeatmapAnnotation(
#           Fenotype = cn,
#           col = list(Fenotype = c("Healthy" = "#C7C7C7", "LumA" = "#FF9E4A", 
#                                   "LumB" = "#67BF5C", "Her2" = "#729ECE", "Basal" = "#ED665D")),
#           annotation_height = 70, show_legend = F, annotation_name_gp = gpar(fontsize = 20)
#         )
# )


