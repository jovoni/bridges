# Robustness of fit(wgd = TRUE) on OV2295 at 500 kb: why are cells placed in
# a sister clone's clade, and does denoising or a different tree method help?
#   Rscript dev/robustness/ov2295_robustness.R [output.rds]
# Runs on the installed package (devtools::install() first).
suppressMessages(library(bridges))
library(data.table)
source("dev/wgd/real_data_prep.R")
source("dev/wgd/wgd_proto.R")   # clone_coherence()

args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args)) args[1] else "ov2295_robustness.rds"
chrs <- as.character(1:22)
rd <- load_ov2295(1)

# ---- Support filter: drop breakpoints carried by fewer than `min_cells` ----
# For each chromosome, a breakpoint is a change of CN between adjacent bins
# in a cell. Breakpoints shared (at the same position) by fewer than
# `min_cells` cells are treated as noise: each cell's profile is re-segmented
# using only its supported breakpoints, and every segment is set to its
# rounded median CN.
support_filter <- function(data, min_cells) {
  dt <- as.data.table(data)
  setorder(dt, chr, start, cell_id)
  out <- list()
  for (ch in unique(dt$chr)) {
    d <- dt[chr == ch]
    M <- as.matrix(dcast(d, cell_id ~ start, value.var = "CN")[, -1])
    cell_ids <- dcast(d, cell_id ~ start, value.var = "CN")$cell_id
    change <- M[, -1, drop = FALSE] != M[, -ncol(M), drop = FALSE]
    supported <- colSums(change) >= min_cells
    for (i in seq_len(nrow(M))) {
      cuts <- which(change[i, ] & supported)
      bounds <- c(0, cuts, ncol(M))
      for (s in seq_len(length(bounds) - 1)) {
        idx <- (bounds[s] + 1):bounds[s + 1]
        M[i, idx] <- round(stats::median(M[i, idx]))
      }
    }
    starts <- as.numeric(colnames(M))
    out[[ch]] <- data.table(cell_id = rep(cell_ids, ncol(M)), chr = ch,
                            start = rep(starts, each = nrow(M)), CN = as.vector(M))
  }
  res <- rbindlist(out)
  res <- merge(res, unique(dt[, .(chr, start, end)]), by = c("chr", "start"))
  dplyr::as_tibble(res[, .(cell_id, chr, start, end, CN)])
}

# ---- Metrics ----
local_misplacement <- function(tree, clones, min_size = 11) {
  cells <- tree$tip.label; n_tip <- length(cells)
  parent <- integer(n_tip + tree$Nnode); parent[tree$edge[, 2]] <- tree$edge[, 1]
  desc <- vector("list", n_tip + tree$Nnode)
  for (i in seq_len(n_tip)) {
    node <- i
    repeat { desc[[node]] <- c(desc[[node]], i); if (node == n_tip + 1) break; node <- parent[node] }
  }
  maj <- vapply(seq_len(n_tip), function(i) {
    node <- i; while (length(desc[[node]]) < min_size) node <- parent[node]
    names(which.max(table(clones[cells[setdiff(desc[[node]], i)]])))
  }, "")
  stats::setNames(maj, cells)
}

evaluate <- function(label, data, k_jitter = 0, tree_func = ape::nj) {
  t <- system.time(res <- fit(data, chromosomes = chrs, alleles = "CN", wgd = TRUE,
                              k_jitter_fix = k_jitter, tree_func = tree_func, verbose = FALSE))
  cpu <- t[["user.self"]] + t[["sys.self"]]
  tree <- res$tree; cells <- tree$tip.label; clones <- rd$clones[cells]
  maj <- local_misplacement(tree, rd$clones)
  mis <- maj != clones
  D <- res$D[cells, cells]; diag(D) <- Inf
  knn <- vapply(seq_along(cells), function(i) mean(clones[order(D[i, ])[1:10]] == clones[i]), 0)
  coh <- clone_coherence(tree, rd$clones)
  sister <- function(a, b) sum(clones[mis] == a & maj[mis] == b) + sum(clones[mis] == b & maj[mis] == a)
  big <- res$wgd$edges[order(res$wgd$edges$n_cells, decreasing = TRUE), ][1:2, ]
  big_desc <- vapply(big$cells, function(x) { t <- table(rd$clones[x]); paste(names(t), t, sep = ":", collapse = " ") }, "")
  data.frame(config = label, misplaced_pct = round(100 * mean(mis), 1),
             CD = sister("C", "D"), AB = sister("A", "B"), GH = sister("G", "H"),
             knn_same_clone = round(mean(knn), 3), clone_f1 = round(weighted.mean(coh$best_f1, coh$n), 3),
             wgd_edges = nrow(res$wgd$edges), wgd_clade_1 = big_desc[1], wgd_clade_2 = big_desc[2],
             cpu_s = round(cpu, 1))
}

results <- list()
add <- function(...) { r <- evaluate(...); print(r[, 1:9], row.names = FALSE); results[[length(results) + 1]] <<- r
                       saveRDS(do.call(rbind, results), out_file) }

RUN_FIRST_ROUND <- FALSE   # jitter 0-2 and support filter: see README
if (RUN_FIRST_ROUND) {
  # 1. k_jitter_fix
  for (k in 0:2) add(sprintf("jitter %d", k), rd$data, k_jitter = k)
  # 2. support filter (alone and with jitter 1)
  for (s in c(3, 5, 10)) {
    filtered <- support_filter(rd$data, s)
    add(sprintf("support >= %d cells", s), filtered)
    add(sprintf("support >= %d cells + jitter 1", s), filtered, k_jitter = 1)
  }
}
# 3. larger jitter
for (k in 3:4) add(sprintf("jitter %d", k), rd$data, k_jitter = k)
# 4. tree methods (bionj refuses distances > 100: rescale; topology unchanged)
bionj_scaled <- function(D) ape::bionj(D / max(D))
add("jitter 2 + bionj", rd$data, k_jitter = 2, tree_func = bionj_scaled)
add("jitter 2 + fastme.bal", rd$data, k_jitter = 2, tree_func = ape::fastme.bal)
add("jitter 2 + fastme.ols", rd$data, k_jitter = 2, tree_func = ape::fastme.ols)

# 5. Reference: plain bin-by-bin L1 distance (no event counting), NJ tree
dt <- as.data.table(rd$data); setorder(dt, cell_id, chr, start)
cells <- unique(dt$cell_id)
M <- matrix(dt$CN, nrow = length(cells), byrow = TRUE, dimnames = list(cells, NULL))
L1 <- as.matrix(stats::dist(M, method = "manhattan"))
clones <- rd$clones[cells]
D <- L1; diag(D) <- Inf
knn <- vapply(seq_along(cells), function(i) mean(clones[order(D[i, ])[1:10]] == clones[i]), 0)
tree <- ape::nj(L1)
maj <- local_misplacement(tree, rd$clones); mis <- maj != clones[tree$tip.label]
coh <- clone_coherence(tree, rd$clones)
cat(sprintf("\nReference, L1 bin distance + NJ (unrooted): misplaced %.1f%% | nearest cells same clone %.3f | clone coherence %.3f\n",
            100 * mean(mis), mean(knn), weighted.mean(coh$best_f1, coh$n)))
