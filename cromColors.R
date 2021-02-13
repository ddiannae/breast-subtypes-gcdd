subtypes.pal <- c("#C7C7C7", "#FF9E4A", "#67BF5C", "#729ECE", "#ED665D")
names(subtypes.pal) <- c("healthy", "luma", "lumb", "her2", "basal")
chrs <- as.character(c(1:22, "X"))
chromosomes.pal <- c("#D909D1", "#0492EE", "#5DA0CB", "#106F35", "#5BD2AE", "#199F41", 
                     "#FE0F43", "#00FFCC", "#F495C5", "#E1BF5D", "#5F166F", "#088ACA",
                     "#41CFE0", "#0F0A71", "#FFFF99", "#B06645", "#800092", "#B925AE",
                     "#B1B719", "#CB97E8", "#130B9E", "#E12B29", "#79A5B9")

names(chromosomes.pal) <- c("22","11","12","13","14","15","16","17","18","19","1" ,"2" ,"3" ,"4" ,"5" ,
                            "6" ,"7" ,"X" ,"8" ,"9" ,"20","10","21")
chromosomes.pal <- chromosomes.pal[chrs]

#cis color #2d20a5, 45,32,165
#trans 255,136,51