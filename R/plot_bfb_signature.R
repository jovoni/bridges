
#' Plot BFB-like Signatures Along the Genome
#'
#' This function visualizes the distribution of BFB-like events (asymmetric CN differences)
#' across genomic segments for a given chromosome and allele. The plot highlights which genomic
#' regions are most frequently affected by potential BFB events in the reconstructed phylogeny.
#'
#' @param res A fitted object returned by the `fit()` function. It must contain the `reconstructions`
#'   field populated via `compute_reconstructions()`.
#' @param chr_of_interest A character string indicating the chromosome to plot (e.g., `"8"` or `"X"`).
#' @param allele_of_interest A character string indicating the allele to analyze, typically `"A"` or `"B"`.
#'
#' @return A `ggplot` object displaying a bar plot where each bar corresponds to a genomic segment,
#'   with its height and color intensity proportional to the number of BFB-like events affecting that segment.
#'
#' @details
#' The function extracts the reconstructed merged profiles and identifies the internal nodes in the
#' phylogeny that show asymmetric copy number changes (i.e., delta > 0). It computes the difference
#' between the left and right child profiles for these nodes, identifies affected segments, and counts
#' the number of times each segment is involved in a BFB-like event. The output is a genomic bar plot
#' showing the number of BFB-like events per segment.
#'
#' @export
plot_bfb_signature = function(res, chr_of_interest, allele_of_interest) {
  rec = res$reconstructions[[chr_of_interest]][[allele_of_interest]]

  deltas = rec$deltas %>% unlist()
  bfb_nodes = names(deltas[deltas > 0])
  bfb_signature = lapply(1:length(bfb_nodes), function(j) {
    rec$merged_profiles[[bfb_nodes[j]]]$left - rec$merged_profiles[[bfb_nodes[j]]]$right
  }) %>% do.call(rbind, .)


  bfb_signature_df = lapply(1:ncol(bfb_signature), function(segment_idx) {
    # Dataframe of Signature
    start_end = unlist(strsplit(colnames(bfb_signature)[segment_idx], ":"))
    start = as.numeric(start_end[1])
    end = as.numeric(start_end[2])
    width = end - start
    df_s = dplyr::tibble(bfb_value = bfb_signature[,segment_idx], node_name = bfb_nodes, segment = segment_idx)

    dplyr::left_join(df_s, dplyr::tibble(segment=segment_idx, start=start, end=end, width=width), by = "segment")
  }) %>% do.call("bind_rows", .)

  bfb_signature_df %>%
    dplyr::mutate(bfb_value = as.numeric(.data$bfb_value != 0)) %>%
    dplyr::group_by(.data$segment) %>%
    dplyr::mutate(bfb_value = sum(.data$bfb_value)) %>%
    dplyr::select(.data$bfb_value, .data$start, .data$end, .data$width) %>%
    dplyr::distinct() %>%
    dplyr::mutate(segment_center = .data$start + .data$width / 2) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$segment_center, y = .data$bfb_value, fill = .data$bfb_value, width = .data$width)) +
    ggplot2::geom_col(position = "identity") +
    ggplot2::scale_fill_gradient(low = "lightblue", high = "darkred") +  # Blue-to-red gradient
    ggplot2::theme_bw() +
    ggplot2::labs(x = paste0("Genomic position (chr", chr_of_interest, ")"), y = "BFB-like events", fill = "BFB count")
}
