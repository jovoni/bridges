
#' Plot Chromosome-wide Heatmap of Copy Number Alterations
#'
#' This function plots a heatmap of copy number alterations (CNA) for a specified chromosome and allele
#' across multiple cells. It supports optional ordering of cells, as well as annotations showing average
#' copy number or gain/loss profiles.
#'
#' @param cna_data A data frame containing copy number data. Must include columns: `cell_id`, chromosome,
#'        bin positions, and a value column specified by `allele`.
#' @param chr Chromosome to plot (e.g., `"1"`, `"X"`).
#' @param allele Column name in `cna_data` containing the copy number values to be plotted: `"CN"`, `"A"` or `"B"`.
#' @param order_heatmap Logical, whether to order the heatmap rows by the number of bins matching a reference value
#'        (2 for total CN, 1 for alleles).
#' @param add_avg_CN_profile Logical, whether to add a barplot annotation showing the average CN profile across cells.
#' @param add_gain_loss_profile Logical, whether to add gain/loss annotation bars (proportion of cells
#'        above/below the baseline: 2 for `"CN"`, 1 for a single allele).
#' @param use_raster Logical, whether to rasterize the heatmap for faster rendering in large datasets.
#' @param raster_quality Integer, quality of the raster image (relevant only if `use_raster = TRUE`).
#'
#' @return A ComplexHeatmap object representing the CNA heatmap, which can be plotted or combined with other heatmaps.
#'
#' @details
#' When \code{order_heatmap = TRUE}, cells are sorted by the number of bins that
#' match the baseline diploid value (2 for total CN, 1 for allele-specific).
#' Cells with the fewest altered bins appear at the top.  The optional annotation
#' bars (\code{add_avg_CN_profile}, \code{add_gain_loss_profile}) are displayed
#' above the heatmap: the average CN profile across cells (bar height = mean CN)
#' and a stacked gain/loss bar (proportion of cells with value above/below
#' baseline at each bin).
#'
#' @examples
#' \dontrun{
#' sim <- bridge_sim(chromosomes = "8", bfb_allele = "8:A",
#'                   max_cells = 128, lambda = 2)
#' # Total CN view, chromosome 8
#' plot_chr_all_heatmap(sim$cna_data, chr = "8", allele = "CN")
#' # Allele A only, no ordering
#' plot_chr_all_heatmap(sim$cna_data, chr = "8", allele = "A",
#'                      order_heatmap = FALSE)
#' }
#'
#' @export
plot_chr_all_heatmap = function(cna_data, chr, allele,
                                order_heatmap = TRUE,
                                add_avg_CN_profile = TRUE,
                                add_gain_loss_profile = TRUE,
                                use_raster = TRUE,
                                raster_quality = 15) {

  cna_chr = cna_data[as.character(cna_data$chr) == as.character(chr), , drop = FALSE]
  if (nrow(cna_chr) == 0) stop("No rows for chromosome '", chr, "' in cna_data.", call. = FALSE)

  ordered_cell_ids = unique(cna_chr$cell_id)
  matrix_data <- prepare_heatmap_matrix(
    data = cna_chr,
    chromosomes_to_plot = chr,
    feature_name = allele,
    ordered_cell_ids = ordered_cell_ids
  )
  mat = matrix_data$matrix

  # Baseline: 2 copies for total CN, 1 for a single allele
  baseline = if (allele == "CN") 2 else 1

  if (order_heatmap) {
    ordered_cell_ids = rowSums(mat == as.character(baseline)) %>% sort() %>% names()
    mat = mat[ordered_cell_ids, , drop = FALSE]
  }

  # Numeric cells x bins matrix (same bin order as `mat`) for the annotations
  numeric_matrix <- tibble_to_matrix(cna_chr, value_column = allele)

  annos = list()
  if (add_avg_CN_profile) {
    annos$AverageCN = ComplexHeatmap::anno_barplot(
      colMeans(numeric_matrix),
      bar_width = 1,
      gp = grid::gpar(fill = "indianred", col = "indianred"),
      border = FALSE
    )
  }

  if (add_gain_loss_profile) {
    annos$Gain = ComplexHeatmap::anno_barplot(
      colMeans(numeric_matrix > baseline),
      bar_width = 1,
      gp = grid::gpar(fill = "indianred", col = "indianred"),
      border = FALSE, ylim = c(0, 1)
    )
    annos$Loss = ComplexHeatmap::anno_barplot(
      -colMeans(numeric_matrix < baseline),
      bar_width = 1,
      gp = grid::gpar(fill = "steelblue", col = "steelblue"),
      border = FALSE, ylim = c(-1, 0)
    )
  }

  top_anno = if (length(annos) > 0) {
    do.call(ComplexHeatmap::HeatmapAnnotation,
            c(annos, list(annotation_name_side = "right")))
  } else {
    NULL
  }

  colvals <- get_colors("CN")

  copynumber_hm = ComplexHeatmap::Heatmap(
    name = "Copy Number",
    as.matrix(mat),
    col = colvals,
    show_column_names = FALSE,
    show_row_names = FALSE,
    use_raster = use_raster,
    raster_quality = raster_quality,
    top_annotation = top_anno,
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )

  copynumber_hm
}

