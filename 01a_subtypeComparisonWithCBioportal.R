library(readr)
library(dplyr)
library(stringr)

our_files <- read_tsv("data/id-filenames-submitterid-caseid-breast-cancer.2019-01-09.tsv")
our_subtypes <- read_tsv("data/id-subtype.tsv")
cbioportal_subtypes <- read_tsv("data/brca_tcga_pan_can_atlas_2018_clinical_data.tsv")
cbioportal_subtypes <- cbioportal_subtypes %>% select(`Patient ID`, Subtype) %>% 
  rename(SubmitterID = `Patient ID`) %>% 
  mutate(Subtype = str_remove(Subtype, "BRCA_"))


our_subtypes <- our_subtypes %>% left_join(our_files, by = "ID") %>% 
  select(SubmitterID, Subtype)
all_subtypes <- our_subtypes %>% 
  inner_join(cbioportal_subtypes, by = "SubmitterID", suffix = c("_our", "_cbio"))
all_subtypes %>% mutate(Equal = Subtype_our == Subtype_cbio) %>% 
  group_by(Equal) %>% tally()
# A tibble: 3 x 2
# Equal     n
# <lgl> <int>
# 1 FALSE   100
# 2 TRUE    554
# 3 NA       68