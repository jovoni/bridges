# Load the OV2295 sample in dev/data as a bridges cna_data table.
# Only clone-labelled cells and autosomes are kept (the synthetic diploid
# reference used for rooting assumes one X and one Y; this sample is female).
# Bins can be merged by `factor` (median state, rounded) to reduce noise.
library(data.table)

load_ov2295 <- function(factor = 1, dir = "dev/data") {
  d <- fread(file.path(dir, "ov2295_cell_cn.csv.gz"),
             select = c("cell_id", "chr", "start", "end", "state"))
  clones <- fread(file.path(dir, "ov2295_clone_clusters.csv.gz"))
  d <- d[cell_id %in% clones$cell_id & chr %in% as.character(1:22)]
  setorder(d, cell_id, chr, start)
  if (factor > 1) {
    d[, bin := (seq_len(.N) - 1L) %/% factor, by = .(cell_id, chr)]
    d <- d[, .(start = min(start), end = max(end), state = round(stats::median(state))),
           by = .(cell_id, chr, bin)]
    d[, bin := NULL]
  }
  d[, CN := as.numeric(state)][, state := NULL]
  list(data = dplyr::as_tibble(d),
       clones = stats::setNames(clones$clone_id, clones$cell_id))
}
