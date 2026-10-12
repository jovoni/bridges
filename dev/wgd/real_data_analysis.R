# Real-data check of the WGD-aware distance on OV2295 (dev/data).
#   Rscript dev/wgd/real_data_analysis.R [output.rds] [figure_dir]
suppressMessages(devtools::load_all(".", quiet = TRUE))
Rcpp::sourceCpp("dev/wgd/wgd_proto.cpp")
source("dev/wgd/wgd_proto.R")
source("dev/wgd/real_data_prep.R")

args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args) >= 1) args[1] else "wgd_real.rds"
fig_dir <- if (length(args) >= 2) args[2] else NULL
chrs <- as.character(1:22)

groups <- list(
  near_diploid = c("A", "B", "C", "D"),
  tetraploid_EF = c("E", "F"),
  triploid_GHI = c("G", "H", "I"),
  all_non_diploid = c("E", "F", "G", "H", "I")
)

group_f1 <- function(tree, clones, members) {
  cells <- names(clones)[clones %in% members]
  parts <- c(lapply(ape::prop.part(tree), function(ix) tree$tip.label[ix]))
  max(vapply(parts, function(p) 2 * length(intersect(p, cells)) / (length(p) + length(cells)), 0))
}

results <- list()
for (f in c(1, 4)) {
  rd <- load_ov2295(f)
  data <- dplyr::bind_rows(rd$data, create_diploid_data(rd$data))
  X <- process_input_data(data, chrs, "CN", 0, 0)
  cls <- estimate_wgd_class(data)
  clones <- rd$clones[names(rd$clones) %in% rownames(X[[1]]$CN)]
  t_comp <- system.time(comp <- wgd_components(X, chrs, "CN"))[["elapsed"]]

  cat(sprintf("\n==== bins x%d (%.1f Mb), %d cells, WGD components %.0fs ====\n",
              f, 0.5 * f, length(clones), t_comp))
  cat("Estimated WGD count per clone:\n")
  print(table(clone = clones, k = cls[names(clones)]))

  for (mode in c("none", "all", "ploidy")) {
    d <- wgd_combine(comp, w = 1, mode = mode, wgd_class = cls)
    tree <- tree_from_distance(d$D)
    coh <- clone_coherence(tree, clones)
    gf1 <- vapply(groups, function(g) group_f1(tree, clones, g), 0)

    kn <- assign_wgd_counts(tree, cls[tree$tip.label])
    rec <- reconstruct_wgd(tree, X, chrs, "CN", kn)
    std <- compute_reconstructions(list(tree = tree, all_input_Xs = X, g_dist_func = "greedy_fast",
                                        b_dist_func = "bfb_fast"), chrs, "CN")
    sets <- wgd_clade_sets(tree, rec)
    comp_tbl <- vapply(sets, function(s) paste(names(table(clones[s])), table(clones[s]), sep = ":", collapse = " "), "")

    nd <- names(clones)[clones %in% groups$near_diploid]
    wg <- names(clones)[!clones %in% groups$near_diploid]
    cat(sprintf("\n-- mode %s --\n", mode))
    cat(sprintf("mean clone best-F1 %.3f | monophyletic clones %d/%d | mean D near-diploid vs others %.1f (within near-diploid %.1f)\n",
                weighted.mean(coh$best_f1, coh$n), sum(coh$monophyletic), nrow(coh),
                mean(d$D[nd, wg]), mean(d$D[nd, nd][upper.tri(d$D[nd, nd])])))
    print(data.frame(clone = coh$clone, n = coh$n, best_f1 = round(coh$best_f1, 3), monophyletic = coh$monophyletic), row.names = FALSE)
    cat("group best-F1:", paste(names(gf1), round(gf1, 3), collapse = " | "), "\n")
    cat(sprintf("WGD edges: %d (root k = %d); events standard %.0f vs WGD-aware %.0f\n",
                length(sets), rec$root_k, sum(std$edge.length), rec$total_events))
    big <- order(lengths(sets), decreasing = TRUE)
    for (i in head(big, 6)) cat(sprintf("  WGD clade of %d cells: %s\n", length(sets[[i]]), comp_tbl[i]))
    if (length(sets) > 6) cat(sprintf("  ... and %d smaller WGD clades (sizes: %s)\n", length(sets) - 6,
                                      paste(sort(lengths(sets[-head(big, 6)]), decreasing = TRUE), collapse = ",")))

    results[[length(results) + 1]] <- list(factor = f, mode = mode, tree = tree, coherence = coh,
                                           group_f1 = gf1, wgd_sets = sets, k_node = kn,
                                           events_standard = sum(std$edge.length),
                                           events_wgd_aware = rec$total_events)

    if (!is.null(fig_dir) && f == 4 && mode %in% c("none", "all")) {
      dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
      ann <- data.frame(cell_id = names(clones), clone = clones,
                        wgd_class = factor(cls[names(clones)]))
      png(file.path(fig_dir, sprintf("ov2295_%s.png", mode)), width = 1800, height = 1200, res = 150)
      ComplexHeatmap::draw(suppressMessages(plot_heatmap(rd$data, tree = tree, chromosomes_to_plot = chrs,
                                                          to_plot = "CN", annotations = ann, branch_length = 1)))
      dev.off()
    }
  }
}
saveRDS(results, out_file)
