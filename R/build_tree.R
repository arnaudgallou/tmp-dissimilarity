####  Init  ####
  {
    library(tidyverse)
    library(V.PhyloMaker2)
  }

####  Build tree  ####
  {
    ptree <- read_csv("data/taxa_df.csv") %>%
      select(species, genus, family) %>%
      phylo.maker(tree = GBOTB.extended.TPL, output.sp.list = FALSE)

    ape::write.tree(ptree$scenario.3, "tree.txt")
  }
