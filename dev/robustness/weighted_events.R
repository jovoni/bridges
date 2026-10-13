# Does a length-saturating event cost make the distance more robust to noise?
# Compares, on OV2295 at 500 kb (no WGD term, which does not affect sister
# clones): current min(G, B), weighted min(G_w, B_w) for several L, and the
# plain bin-by-bin L1 distance.
#   Rscript dev/robustness/weighted_events.R
suppressMessages(library(bridges))
library(data.table)
Rcpp::sourceCpp("dev/robustness/weighted_events.cpp")
source("dev/wgd/real_data_prep.R")
source("dev/wgd/wgd_proto.R")
src <- readLines("dev/robustness/ov2295_robustness.R")
eval(parse(text = src[grep("^local_misplacement <- function", src):(grep("^evaluate <- function", src) - 1)]))

chrs <- as.character(1:22)
rd <- load_ov2295(1)
bin <- 5e5

score <- function(label, D) {
  cells <- rownames(D); clones <- rd$clones[cells]
  Dk <- D; diag(Dk) <- Inf
  knn <- vapply(seq_along(cells), function(i) mean(clones[order(Dk[i, ])[1:10]] == clones[i]), 0)
  tree <- ape::nj(D)
  maj <- local_misplacement(tree, rd$clones); mis <- maj != clones[tree$tip.label]
  coh <- clone_coherence(tree, rd$clones)
  sister <- function(a, b) sum(clones[tree$tip.label][mis] %in% c(a, b) & maj[mis] %in% c(a, b))
  cat(sprintf("%-34s misplaced %5.1f%% | C/D %3d  A/B %3d  G/H %3d | nearest cells same clone %.3f | clone coherence %.3f\n",
              label, 100 * mean(mis), sister("C", "D"), sister("A", "B"), sister("G", "H"),
              mean(knn), weighted.mean(coh$best_f1, coh$n)))
}

for (jit in c(0, 2)) {
  X <- bridges:::process_input_data(rd$data, chrs, "CN", jit, 0)
  cells <- rownames(X[[1]]$CN)
  for (L in c(0, 2, 5, 10, 20)) {
    D <- matrix(0, length(cells), length(cells), dimnames = list(cells, cells))
    for (ch in chrs) {
      M <- X[[ch]]$CN[cells, , drop = FALSE]
      se <- do.call(rbind, strsplit(colnames(M), ":"))
      widths <- (as.numeric(se[, 2]) - as.numeric(se[, 1]) + 1) / bin
      D <- D + weighted_min_distance(M, widths, L)
    }
    score(sprintf("jitter %d, event cost min(1, len/%d)", jit, L), D)
  }
}
dt <- as.data.table(rd$data); setorder(dt, cell_id, chr, start)
cells <- unique(dt$cell_id)
M <- matrix(dt$CN, nrow = length(cells), byrow = TRUE, dimnames = list(cells, NULL))
score("reference: L1 bin distance", as.matrix(stats::dist(M, method = "manhattan")))
