#' Centromere positions and chromosome arm boundaries (hg19/hg38)
#'
#' Real centromere positions, precomputed p-/q-arm boundaries, and total
#' chromosome lengths for the 24 standard human chromosomes (1-22, X, Y)
#' under both hg19 (GRCh37) and hg38 (GRCh38). Used internally by
#' \code{\link{split_chr_arms}} (via the \code{genome} argument) and by
#' \code{\link{bridge_sim}} to size simulated chromosomes.
#'
#' @format A data frame with 48 rows (24 chromosomes x 2 genome builds) and
#'   columns:
#' \describe{
#'   \item{genome}{\code{"hg19"} or \code{"hg38"}.}
#'   \item{chr}{Bare chromosome id, \code{"1".."22"}, \code{"X"}, \code{"Y"}
#'     (no \code{"chr"} prefix).}
#'   \item{chr_length}{Total chromosome length in bp.}
#'   \item{centromere_start,centromere_end}{1-based bounds of the centromeric
#'     region.}
#'   \item{centromere_mid}{Midpoint of the centromeric region -- the value
#'     \code{split_chr_arms()} thresholds on.}
#'   \item{p_arm_start,p_arm_end,q_arm_start,q_arm_end}{Arm boundaries.}
#' }
#' @source UCSC Genome Browser \code{gap}/\code{centromeres} and
#'   \code{chrom.sizes} tracks for hg19 and hg38, downloaded 2026-07-22. See
#'   \code{data-raw/centromeres.R} for the exact source URLs and transform.
#' @examples
#' centromeres[centromeres$genome == "hg38" & centromeres$chr == "1", ]
"centromeres"

utils::globalVariables("centromeres")

.centromere_vector_for <- function(genome) {
  genome <- match.arg(genome, choices = c("hg19", "hg38"))
  tbl <- centromeres[centromeres$genome == genome, ]
  stats::setNames(tbl$centromere_mid, tbl$chr)
}

.chr_lengths_for <- function(genome, chromosomes) {
  genome <- match.arg(genome, choices = c("hg19", "hg38"))
  tbl <- centromeres[centromeres$genome == genome, ]
  lens <- stats::setNames(tbl$chr_length, tbl$chr)
  missing <- setdiff(as.character(chromosomes), names(lens))
  if (length(missing)) {
    stop("No chr_length available for chromosome(s): ", paste(missing, collapse = ", "),
         " under genome = '", genome, "'. bridges bundles data for chromosomes ",
         "1-22, X, Y only.", call. = FALSE)
  }
  lens[as.character(chromosomes)]
}
