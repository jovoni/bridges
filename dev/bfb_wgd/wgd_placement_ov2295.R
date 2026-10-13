# Where are the WGD events placed on the OV2295 tree (500 kb, jitter 2), and
# what are the ancestral profiles just before and just after each WGD?
#   Rscript dev/bfb_wgd/wgd_placement_ov2295.R [figure_dir]
suppressMessages(library(bridges))
source("dev/wgd/real_data_prep.R")
args <- commandArgs(trailingOnly = TRUE)
fig_dir <- if (length(args)) args[1] else "dev/data/figures"
chrs <- as.character(1:22); rd <- load_ov2295(1)
res <- fit(rd$data, chromosomes = chrs, alleles = "CN", k_jitter_fix = 2, wgd = TRUE, verbose = FALSE)
tree <- res$tree; n_tip <- length(tree$tip.label); clones <- rd$clones
tips_of <- function(node) if (node <= n_tip) tree$tip.label[node] else ape::extract.clade(tree, node)$tip.label
comp <- function(cells) { t <- table(clones[cells]); paste(names(t), t, sep = ":", collapse = " ") }
G <- function(a, b) bridges:::greedy_distance_cpp(as.numeric(a), as.numeric(b), 2)$cost

# Profile of a node for one chromosome (tips: observed; internal: reconstructed)
node_profile <- function(chr, node) {
  r <- res$reconstructions[[chr]][["CN"]]
  if (node <= n_tip) res$all_input_Xs[[chr]][["CN"]][tree$tip.label[node], ] else r$internal_nodes[[as.character(node)]]
}

e <- res$wgd$edges[order(res$wgd$edges$n_cells, decreasing = TRUE), ]
cat(sprintf("WGD counts per cell: %s | WGD above the root: %d\n",
            paste(names(table(res$wgd$counts)), table(res$wgd$counts), sep = ":", collapse = " "), res$wgd$root_count))
cat(sprintf("%d WGD branches:\n", nrow(e)))
for (i in seq_len(nrow(e))) {
  sib <- setdiff(tree$edge[tree$edge[, 1] == e$parent[i], 2], e$child[i])
  cat(sprintf("  %3d cells (%s) | sister clade: %d cells (%s)\n", e$n_cells[i], comp(e$cells[[i]]),
              length(unlist(lapply(sib, tips_of))), comp(unlist(lapply(sib, tips_of)))))
}

profiles <- list()
for (i in 1:2) {
  parent <- e$parent[i]; child <- e$child[i]; cells <- e$cells[[i]]
  label <- if (any(clones[cells] %in% c("G", "H", "I"))) "G/H/I" else "E/F"
  cat(sprintf("\n== WGD on the %s branch (%d cells) ==\n", label, length(cells)))
  cat("events per chromosome: diploid -> pre-WGD ancestor | 2 x pre-WGD -> post-WGD ancestor (on the WGD branch) | post-WGD ancestor -> clade consensus\n")
  rows <- list()
  for (chr in chrs) {
    pre <- node_profile(chr, parent); post <- node_profile(chr, child)
    X <- res$all_input_Xs[[chr]][["CN"]][cells, , drop = FALSE]
    cons <- round(apply(X, 2, stats::median))
    seg <- do.call(rbind, strsplit(colnames(res$all_input_Xs[[chr]][["CN"]]), ":"))
    profiles[[length(profiles) + 1]] <- data.frame(wgd = label, chr = chr, start = as.numeric(seg[, 1]), end = as.numeric(seg[, 2]),
                                                   pre_wgd = pre, doubled_pre = 2 * pre, post_wgd = post, clade_consensus = cons)
    rows[[chr]] <- c(before = G(rep(2, length(pre)), pre), on_wgd_branch = G(2 * pre, post), within_clade = G(post, cons))
  }
  tab <- do.call(rbind, rows)
  print(t(tab))
  cat(sprintf("totals: before WGD %d | on WGD branch after doubling %d | within clade %d\n",
              sum(tab[, 1]), sum(tab[, 2]), sum(tab[, 3])))
}
prof <- do.call(rbind, profiles)
saveRDS(prof, file.path(fig_dir, "ov2295_wgd_ancestral_profiles.rds"))

# Chromosome 6 close-up for the G/H/I WGD
p6 <- prof[prof$wgd == "G/H/I" & prof$chr == "6", ]
cells <- e$cells[[1]]
X6 <- res$all_input_Xs[["6"]][["CN"]][cells, , drop = FALSE]
p6$frac_cells_cn_ge_6 <- round(colMeans(X6 >= 6), 2)
cat("\nchromosome 6, G/H/I WGD (segments, Mb):\n")
print(transform(p6[, c("start", "end", "pre_wgd", "doubled_pre", "post_wgd", "clade_consensus", "frac_cells_cn_ge_6")],
                start = round(start / 1e6, 1), end = round(end / 1e6, 1)), row.names = FALSE)

# Figure: genome-wide profiles before / after each WGD
library(ggplot2)
long <- rbind(
  data.frame(prof[, c("wgd", "chr", "start", "end")], profile = "2 x pre-WGD ancestor", CN = prof$doubled_pre),
  data.frame(prof[, c("wgd", "chr", "start", "end")], profile = "post-WGD ancestor", CN = prof$post_wgd),
  data.frame(prof[, c("wgd", "chr", "start", "end")], profile = "clade consensus", CN = prof$clade_consensus))
long$chr <- factor(long$chr, levels = chrs)
g <- ggplot(long, aes(x = start / 1e6, xend = end / 1e6, y = CN, yend = CN, colour = profile)) +
  geom_segment(linewidth = 1, position = position_nudge(y = 0)) +
  facet_grid(wgd ~ chr, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = c("2 x pre-WGD ancestor" = "grey50", "post-WGD ancestor" = "steelblue",
                                 "clade consensus" = "firebrick")) +
  scale_x_continuous(breaks = NULL) +
  labs(x = NULL, y = "Copy number", colour = NULL,
       title = "OV2295 (500 kb, jitter 2): reconstructed profiles around each WGD") +
  theme_bw() + theme(legend.position = "bottom", panel.spacing = unit(1, "pt"),
                     strip.text.x = element_text(size = 7))
ggsave(file.path(fig_dir, "ov2295_wgd_ancestral_profiles.png"), g, width = 16, height = 6, dpi = 150)
cat("\nsaved", file.path(fig_dir, "ov2295_wgd_ancestral_profiles.png"), "\n")
