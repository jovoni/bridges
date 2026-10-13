# Event model vs the current greedy/BFB model, both with WGD terms.
# Same scenarios and seeds as run_experiments_hard.R, noise "none" and "low".
#   Rscript dev/wgd/run_event_model.R [output.rds]
src <- readLines("dev/wgd/run_experiments_hard.R")
eval(parse(text = src[1:(grep("# ---- Evaluation", src, fixed = TRUE) - 1)]))
Rcpp::sourceCpp("dev/wgd/event_model.cpp")
source("dev/wgd/event_model.R")
args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args)) args[1] else "wgd_event_model.rds"

methods <- list(
  greedy_no_wgd = list(model = "greedy", mode = "none"),
  greedy_wgd    = list(model = "greedy", mode = "all"),
  event_wgd     = list(model = "event", mode = "all", use_bfb = FALSE),
  event_bfb_wgd = list(model = "event", mode = "all", use_bfb = TRUE)
)

results <- list()
for (sc in names(scenarios)) for (rep in seq_along(scenarios[[sc]])) {
  s <- scenarios[[sc]][[rep]]
  for (nl in c("none", "low")) {
    set.seed(s$seed + 7)
    noisy <- add_noise(s$data, noise_levels[[nl]][1], noise_levels[[nl]][2])
    for (alleles in list(c("A", "B"), "CN")) {
      data <- dplyr::bind_rows(noisy, create_diploid_data(noisy))
      X <- process_input_data(data, chrs, alleles, 0, 0)
      cls <- estimate_wgd_class(data)
      cells <- names(s$n_wgd); k <- s$n_wgd[cells]
      true_sets <- c(s$events, if (any(k == 2)) list(k2 = cells[k == 2]))
      t_g <- system.time(comp_g <- wgd_components(X, chrs, alleles))[["elapsed"]]
      t_e <- system.time(comp_e <- event_components(X, chrs, alleles, use_bfb = FALSE))[["elapsed"]]
      t_eb <- system.time(comp_eb <- event_components(X, chrs, alleles, use_bfb = TRUE))[["elapsed"]]
      for (m in names(methods)) {
        cfg <- methods[[m]]
        comp <- if (cfg$model == "greedy") comp_g else if (isTRUE(cfg$use_bfb)) comp_eb else comp_e
        tree <- tree_from_distance(wgd_combine(comp, 1, cfg$mode, cls)$D)
        kn <- assign_wgd_counts(tree, cls[tree$tip.label])
        rec <- if (cfg$model == "greedy") reconstruct_wgd(tree, X, chrs, alleles, kn) else
          reconstruct_wgd_event(tree, X, chrs, alleles, kn)
        sc_lab <- score_wgd_labels(true_sets, wgd_clade_sets(tree, rec))
        results[[length(results) + 1]] <- data.frame(
          scenario = sc, rep = rep, noise = nl, data = paste(alleles, collapse = "+"), method = m,
          rf = phangorn::RF.dist(s$tree, tree, normalize = TRUE),
          label_jaccard = if (length(true_sets)) mean(sc_lab$best_jaccard) else NA,
          n_true_wgd = sc_lab$n_true, n_pred_wgd = sc_lab$n_pred,
          distance_seconds = if (cfg$model == "greedy") t_g else if (isTRUE(cfg$use_bfb)) t_eb else t_e)
      }
    }
  }
  message(sc, " rep ", rep, " done")
  saveRDS(do.call(rbind, results), out_file)
}
