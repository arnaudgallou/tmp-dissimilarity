mtx_filter <- function(x, y) {
  x[, colnames(x) %in% intersect(colnames(x), y)]
}

first_word <- function(x) {
  string_extract(x, r"{^[^\pL]*\K[\pL]+}")
}

string_extract <- function(string, pattern) {
  x <- regexpr(pattern, string, perl = TRUE)
  out <- rep(NA, length(string))
  out[x != -1] <- regmatches(string, x)
  out
}

tree_to_genus <- function(tree) {
  tree$tip.label <- first_word(tree$tip.label)
  keep.tip(tree, unique(tree$tip.label))
}

prune <- function(tree, taxa) {
  UseMethod("prune")
}

prune.phylo <- function(tree, taxa) {
  keep.tip(tree, intersect(tree$tip.label, taxa))
}

extend <- function(x, n = 2) {
  map_dfr(1:n, ~ mutate(x, grids = paste0("site_", .x)))
}

sample_by_site <- function(x, n) {
  x <- group_by(x, site)
  x <- slice_sample(x, n = n)
  ungroup(x)
}

semi_difference <- function(x, n) {
  x <- group_split(x, grids)
  imap_dfr(x, ~ {
    i <- sample.int(nrow(.x), n)
    slice(.x, -i)
  })
}

difference <- function(x, n, level = c("species", "genus")) {
  level <- match.arg(level)
  taxa <- sample(x[[level]], n)
  x <- group_split(x, grids)
  imap_dfr(x, ~ {
    if (.y == 1) {
      filter(.x, !.data[[level]] %in% taxa)
    } else {
      filter(.x, .data[[level]] %in% taxa)
    }
  })
}

calc_phylo_dist <- function(data, tree) {
  common_sp <- intersect(tree$tip.label, colnames(data))
  tree <- keep.tip(tree, common_sp)
  x <- mtx_filter(data, common_sp)
  x <- phylobeta(x, tree, index.family = "sorensen")
  c(x$phylo.beta.sor)
}

get_taxa_tbl <- function(x, .genus, .species, level = c("species", "genus")) {
  level <- match.arg(level)
  if (!missing(.species) && level == "species") {
    if (length(.genus) > 1) {
      .genus <- paste(.genus, collapse = "|")
    }
    pattern <- sprintf("(?i)(?:%s)_(?:%s)", .genus, .species)
    .expr <- expr(str_detect(.data[["species"]], pattern))
  } else {
    .expr <- expr(tolower(genus) %in% tolower(.genus))
    if (level == "genus") {
      x <- distinct(x, family, genus)
      x <- mutate(x, species = genus)
    }
  }
  filter(x, !!.expr)
}

get_taxa_list <- function(x, level = c("genus", "species"), within) {
  level <- match.arg(level)
  if (level == "genus") {
    return(unique(x$genus))
  }
  if (missing(within)) {
    return(x$species)
  }
  x <- filter(x, tolower(genus) %in% tolower(within))
  x$species
}
