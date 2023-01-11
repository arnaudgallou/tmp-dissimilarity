####  Init  ####
  {
    library(tidyverse)
    library(ape)
    source("R/helpers.R")

    mdf <- read_csv("data/taxa_df.csv")
  }

####  Helpers  ####
  {
    prune <- function(x, y) {
      UseMethod("prune")
    }

    prune.phylo <- function(x, y) {
      keep.tip(x, y)
    }

    prune.dgCMatrix <- function(x, y) {
      x[, colnames(x) %in% y]
    }

    drop_absent <- function(x) {
      UseMethod("drop_absent")
    }

    drop_absent.dgCMatrix <- function(x) {
      x[, Matrix::colSums(x) > 0]
    }

    make_distance_df <- function(data, tree, pairs, id, drop_absent = TRUE) {
      apply(pairs, 2, function(i) {
        list(
          item = i[1],
          item_2 = i[2],
          distance = calc_mntd(data[i, ], tree, drop_absent)
        )
      })
    }

    equalize <- function(x, tree, drop_absent) {
      if (drop_absent) {
        x <- drop_absent(x)
      }
      # cat("Community matrix:\n")
      # print(x)
      # cat("\n")
      tip_labels <- tree$tip.label
      taxa <- tip_labels[tip_labels %in% colnames(x)]
      list(
        matrix = prune(x, taxa),
        tree = prune(tree, taxa)
      )
    }

    get_communities <- function(x) {
      lapply(seq(nrow(x)), function (i) {
        colnames(x)[x[i, ] == 1]
      })
    }

    mntd <- function(x) {
      # cat("Distance matrix:\n")
      # print(x)
      # cat("\n")
      out <- lapply(1:2, function(i) {
        apply(x, i, min)
      })
      mean(unlist(out))
    }

    calc_mntd <- function(x, tree, drop_absent) {
      data <- equalize(x, tree, drop_absent)
      communities <- get_communities(data$matrix)
      distances <- ape::cophenetic.phylo(data$tree)
      out <- distances[communities[[1]], communities[[2]], drop = FALSE]
      if (length(out) == 1L) {
        return(as.numeric(out))
      }
      mntd(out)
    }
  }

#### Tree ####
  {
    ptree <- read.tree("data/tree.txt")
  }

####  General case  ####
  {
    # to list available genera
    get_taxa_list(mdf)

    set.seed("011023")

    communities <- mdf %>%
      get_taxa_tbl(c("erigeron", "rosa")) %>%
      extend(n = 4) %>%
      semi_difference(n = 300) %>%
      to_sparse_matrix(row_names = "grids", col_names = "species")

    pairs <- combn(rownames(communities), m = 2)

    tibble(data = make_distance_df(
      communities,
      tree = ptree,
      pairs = pairs,
      drop_absent = TRUE
    )) %>%
      unnest_wider(col = data)
  }

####  Influence of the number of taxa in matrix on MNTD  ####
  {
    communities <- mdf %>%
      get_taxa_tbl("erigeron", "vagus|compositus|clokeyi|coulteri|mancus") %>%
      extend(n = 4)

    set.seed("011023")

    df <- semi_difference(communities, n = 3)

    smat <- to_sparse_matrix(df, row_names = "grids", col_names = "species")
    smat <- smat[c(1, 3), ]

    pairs <- combn(rownames(smat), m = 2)

    tibble(data = make_distance_df(
      smat,
      tree = ptree,
      pairs = pairs,
      drop_absent = TRUE
    )) %>%
      unnest_wider(col = data)
  }
