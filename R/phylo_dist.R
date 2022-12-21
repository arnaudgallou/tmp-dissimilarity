####  Init  ####
  {
    library(tidyverse)
    library(phyloregion)
    library(ape)
    source("R/helpers.R")
  }

ptree <- read.tree("data/tree.txt")

taxa <- c(
  c(
    get_taxa_list(mdf, "species", within = "erigeron")[c(1, 7)],
    get_taxa_list(mdf, "species", within = "rosa")[3]
  ),
  c(
    get_taxa_list(mdf, "species", within = "erigeron")[2],
    get_taxa_list(mdf, "species", within = "rosa")[3]
  )
)

# matrix of actual phylogenetic distances between taxa
phylo_dist <- cophenetic(prune(ptree, taxa))

mtx <- tibble(
  grids = rep(paste0("site", 1:2), times = c(3, 2)),
  species = taxa
) %>%
  long2sparse()

mean(phylo_dist[
  colnames(mtx)[mtx["site1", ] == 1],
  colnames(mtx)[mtx["site2", ] == 1]
])
