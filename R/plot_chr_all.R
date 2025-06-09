
#' Plot Chromosome-wide Heatmap of Copy Number Alterations
#'
#' This function plots a heatmap of copy number alterations (CNA) for a specified chromosome and allele
#' across multiple cells. It supports optional ordering of cells, as well as annotations showing average
#' copy number or gain/loss profiles.
#'
#' @param cna_data A data frame containing copy number data. Must include columns: `cell_id`, chromosome,
#'        bin positions, and a value column specified by `allele`.
#' @param chr Chromosome to plot (e.g., `"1"`, `"X"`).
#' @param allele Column name in `cna_data` containing the copy number values to be plotted (e.g., `"CN"`, `"cn_a"`).
#' @param order_heatmap Logical, whether to order the heatmap rows by the number of bins matching a reference value
#'        (2 for total CN, 1 for alleles).
#' @param add_avg_CN_profile Logical, whether to add a barplot annotation showing the average CN profile across cells.
#' @param add_gain_loss_profile Logical, whether to add gain/loss annotation bars (proportion of cells with gain/loss).
#' @param use_raster Logical, whether to rasterize the heatmap for faster rendering in large datasets.
#' @param raster_quality Integer, quality of the raster image (relevant only if `use_raster = TRUE`).
#'
#' @return A ComplexHeatmap object representing the CNA heatmap, which can be plotted or combined with other heatmaps.
#'
#' @export
plot_chr_all_heatmap = function(cna_data, chr, allele,
                                order_heatmap = TRUE,
                                add_avg_CN_profile = TRUE,
                                add_gain_loss_profile = TRUE,
                                use_raster = TRUE,
                                raster_quality = 15) {

  ordered_cell_ids = unique(cna_data$cell_id)
  matrix_data <- prepare_heatmap_matrix(
    data = cna_data,
    chromosomes_to_plot = chr,
    feature_name = allele,
    ordered_cell_ids = ordered_cell_ids
  )
  mat = matrix_data$matrix

  if (order_heatmap) {
    target_val = ifelse(allele == "CN", "2", "1")
    ordered_cell_ids = rowSums(mat == target_val) %>% sort() %>% names()
    mat = mat[ordered_cell_ids, ]
  }

  if (add_avg_CN_profile) {
    numeric_matrix <- tibble_to_matrix(cna_data, chr, value_column = allele)
    avg_CN_profile = colMeans(numeric_matrix)
    top_anno <- ComplexHeatmap::HeatmapAnnotation(
      AverageCN = ComplexHeatmap::anno_barplot(
        avg_CN_profile,
        bar_width = 1,
        beside = FALSE,
        gp = grid::gpar(fill = "indianred", col="indianred"),  # Colors for below and above x
        border = FALSE
      ),
      annotation_name_side = "right"
    )
  }

  if (add_gain_loss_profile) {
    gain_prop <- colMeans(mat > 1) # Proportion of amplified bins
    loss_prop <- colMeans(mat < 1) # Proportion of lost bins
    # Create a stacked bar annotation
    top_anno <- ComplexHeatmap::HeatmapAnnotation(
      Gain = ComplexHeatmap::anno_barplot(
        gain_prop,
        bar_width = 1,
        beside = FALSE,
        gp = grid::gpar(fill = "indianred", col="indianred"),  # Colors for below and above x
        border = FALSE, ylim = c(0,1)
      ),
      Loss = ComplexHeatmap::anno_barplot(
        -loss_prop,
        beside = FALSE,
        bar_width = 1,
        gp = grid::gpar(fill = "steelblue", col="steelblue"),  # Colors for below and above x
        border = FALSE, ylim = c(-1,0)
      ),
      annotation_name_side = "right"
    )
  }

  colvals <- get_colors("CN")

  if (add_gain_loss_profile | add_avg_CN_profile) {
    copynumber_hm = ComplexHeatmap::Heatmap(
      name = "Copy Number",
      as.matrix(mat),
      col = colvals,
      show_column_names = FALSE,
      show_row_names = FALSE,
      use_raster = use_raster,
      raster_quality = raster_quality,
      top_annotation = top_anno,
      cluster_rows = F,
      cluster_columns = F
    )
  } else {
    copynumber_hm = ComplexHeatmap::Heatmap(
      name = "Copy Number",
      as.matrix(mat),
      col = colvals,
      show_column_names = FALSE,
      show_row_names = FALSE,
      use_raster = use_raster,
      raster_quality = raster_quality,
      cluster_rows = F,
      cluster_columns = F
    )
  }

  copynumber_hm
}

