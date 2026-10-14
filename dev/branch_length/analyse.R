# Summaries and figure for events_vs_depth.R.
#   Rscript dev/branch_length/analyse.R events_vs_depth.rds [figure.png]
suppressMessages({ library(dplyr); library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
d <- readRDS(args[1])
fig <- if (length(args) > 1) args[2] else NULL

# Per simulation and method: correlation and scale between truth and estimates
per_run <- d %>%
  group_by(scenario, lambda, rep, fit_wgd) %>%
  summarise(n = n(), true_mean = mean(true_total),
            r_depth = cor(true_total, depth), rho_depth = cor(true_total, depth, method = "spearman"),
            ratio_depth = mean(depth) / mean(true_total),
            r_D = cor(true_total, D_dipl), ratio_D = mean(D_dipl) / mean(true_total),
            .groups = "drop")
cat("== Mean over replicates (r = Pearson, rho = Spearman, ratio = mean estimate / mean truth) ==\n")
print(as.data.frame(per_run %>% group_by(scenario, lambda, fit_wgd) %>%
  summarise(true_events = round(mean(true_mean), 1),
            r_depth = round(mean(r_depth), 3), rho_depth = round(mean(rho_depth), 3),
            ratio_depth = round(mean(ratio_depth), 2),
            r_D = round(mean(r_D), 3), ratio_D = round(mean(ratio_D), 2), .groups = "drop")))

# How much each true event type adds to the depth: depth ~ amp + del + bfb + wgd,
# within each scenario/lambda/method (replicate fixed effects)
cat("\n== Depth added per true event, by type (lm with replicate intercepts) ==\n")
coefs <- d %>% group_by(scenario, lambda, fit_wgd) %>% group_modify(function(x, ...) {
  keep <- c("amp", "del", "bfb", "wgd")[vapply(c("amp", "del", "bfb", "wgd"), function(t) var(x[[t]]) > 0, TRUE)]
  m <- lm(reformulate(c(keep, "factor(rep)"), "depth"), data = x)
  as.data.frame(t(round(coef(m)[keep], 2)))
})
print(as.data.frame(coefs))

if (!is.null(fig)) {
  x <- filter(d, lambda == 1) %>% mutate(method = paste0("fit(wgd = ", fit_wgd, ")"))
  p <- ggplot(x, aes(true_total, depth)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_jitter(aes(colour = factor(wgd)), width = 0.2, height = 0, alpha = 0.35, size = 0.8) +
    facet_grid(method ~ scenario) +
    labs(x = "true events in the cell's lineage", y = "root-to-tip length of the fitted tree",
         colour = "true WGDs", title = "lambda = 1, 5 replicates per panel; dashed: y = x") +
    theme_bw()
  ggsave(fig, p, width = 11, height = 6, dpi = 150)
}
