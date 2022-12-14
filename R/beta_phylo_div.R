####  Init  ####
  {
    library(tidyverse)
    library(phyloregion)
    library(ape)
    source("R/helpers.R")
  }

####  Data  ####
  {
    mdf <- read_csv("data/taxa_df.csv")
  }

####  Dissimilarity  ####
  {
    get_taxa_list(mdf, "genus")
    get_taxa_list(mdf, "species", within = "erigeron")

    # · Species level ----
      {
        ptree <- read.tree("data/tree.txt")

        communities <- mdf %>%
          get_taxa_tbl("erigeron", "vagus|compositus|clokeyi|coulteri|mancus") %>%
          extend()

        communities %>%
          # `difference` makes communities different (no shared taxa)
          # `n` determines the number of unshared taxa between the 2 communities
          difference(n = 2) %>%
          long2sparse() %>% # remove the pipe to see community structures
          calc_phylo_dist(ptree)

        communities %>%
          # `semi_difference` makes communities partially similar (some taxa are
          # found in both communities while others are not)
          # `n` determines the number of taxa unique to one community
          semi_difference(n = 1) %>%
          long2sparse() %>%
          calc_phylo_dist(ptree)
      }

    # · Genus level ----
      {
        ptree <- tree_to_genus(ptree)
        plot(ptree)

        communities <- mdf %>%
          get_taxa_tbl(c("spiraea", "cydonia"), level = "genus") %>%
          extend()

        communities %>%
          difference(n = 1, level = "genus") %>%
          long2sparse() %>%
          calc_phylo_dist(ptree)

        communities %>%
          semi_difference(n = 1) %>%
          long2sparse() %>%
          calc_phylo_dist(ptree)
      }
  }
