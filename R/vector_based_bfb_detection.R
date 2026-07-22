# ============================================================
# bfb_detect_batch.R
#
# Glue code between `bridges` simulation output and `bfbtools`
# detection. Deliberately kept OUTSIDE both packages: bfbtools stays a
# clean, sim-agnostic BFB detection library; bridges stays a clean
# simulator. This file is the (small, easy-to-audit) bridge between the
# two -- drop it into your analysis scripts, or fold it into `bridges`
# itself as a detection helper, whichever fits your project better.
#
# Three functions:
#   - bfb_detect_vectors(vec_list, ...)  : core batch logic over a plain
#     list of count-vectors. Sim-agnostic, fully tested below against
#     bfbtools directly.
#   - split_chr_arms(cna_data, ...)       : adds a p/q arm label to a
#     cna_data tibble, needed because ascending genomic coordinate is
#     telomere-to-centromere on the p-arm but centromere-to-telomere on
#     the q-arm -- only one of the two needs reversing before BFB
#     count-vector extraction.
#   - bfb_detect_batch(cna_data, ...)     : the bridges-facing wrapper,
#     extracting count-vectors from a cna_data tibble (e.g.
#     sim$cna_data from a bridge_sim() result) and calling
#     bfb_detect_vectors(). The chr/allele-splitting and dedup logic is
#     tested below end-to-end against a minimal stand-in for
#     tibble_to_matrix() with the assumed interface -- but
#     NOT against the real bridges implementation (not available in
#     this environment). See the "Assumptions" note on bfb_detect_batch
#     and please sanity-check column names against your actual
#     bridge_sim() output.
# ============================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Run BFB detection over a list of count-vectors
#'
#' Core batch logic, decoupled from any particular simulator: given a
#' list of count-vectors (one per cell/allele/whatever grouping you
#' have), runs \code{bfbtools::admits_bfb()} on each, and
#' \code{bfbtools::nearest_bfb()} for the ones that don't admit BFB
#' exactly -- exploiting exact-duplicate profiles across entries (very
#' common for related cells in a simulated lineage) so each distinct
#' count-vector is only ever computed once, however many entries share
#' it.
#'
#' TRIVIAL-VECTOR FILTER: "does this count-vector admit a BFB schedule"
#' is technically true even for a schedule with *zero* duplication
#' cycles -- the untouched baseline chromosome \code{x0 = (1,2,...,k)} is
#' itself a (trivial) BFB schedule. So any count-vector where every
#' segment sits at copy number 1 (e.g. a flat region interrupted by a
#' plain deletion, with no amplification anywhere) will "admit BFB"
#' without representing anything BFB-like biologically -- BFB is
#' fundamentally an amplification mechanism, and every real cycle is an
#' inverted *duplication*, which necessarily pushes at least one
#' segment above copy number 1. \code{min_max_count} (default 2)
#' enforces that precondition: vectors whose maximum count doesn't
#' exceed it are treated as uninformative (NA results), the same
#' bucket as length <= 1 or NA-containing vectors. Set to 1 to disable
#' this filter and restore the raw (weaker) admits-BFB semantics.
#'
#' @param vec_list a list of integer count-vectors. NULL, length <= 1,
#'   NA-containing, or (by default) max-count-1 entries are treated as
#'   uninformative (NA results).
#' @param model,min_weight passed to \code{bfbtools::nearest_bfb()}.
#' @param min_max_count minimum peak copy number required for a vector
#'   to be considered informative (default 2; see Details above).
#' @return a tibble with one row per element of vec_list: n_segments,
#'   max_count, admits_bfb, bfb_distance (= -log(weight); 0 if
#'   admits_bfb is TRUE), nearest_weight, and n_unique (diagnostic: how
#'   many distinct count-vectors were actually computed).
#' @export
bfb_detect_vectors <- function(vec_list, model = c("poisson", "none"), min_weight = 1e-6,
                                min_max_count = 2) {
  model <- match.arg(model)

  is_trivial <- function(v) {
    is.null(v) || length(v) <= 1 || anyNA(v) || max(v) < min_max_count
  }

  keys <- vapply(vec_list, function(v) {
    if (is_trivial(v)) return(NA_character_)
    paste(v, collapse = ",")
  }, character(1))

  unique_keys <- unique(stats::na.omit(keys))
  cache <- new.env(parent = emptyenv())

  for (key in unique_keys) {
    v <- as.integer(strsplit(key, ",")[[1]])
    entry <- tryCatch({
      admits <- bfbtools::admits_bfb(v)
      if (isTRUE(admits)) {
        list(admits_bfb = TRUE, bfb_distance = 0, nearest_weight = 1)
      } else {
        nb <- suppressWarnings(bfbtools::nearest_bfb(v, model = model, min_weight = min_weight))
        if (is.null(nb) || (length(nb) == 1 && is.na(nb))) {
          list(admits_bfb = FALSE, bfb_distance = NA_real_, nearest_weight = NA_real_)
        } else {
          list(admits_bfb = FALSE, bfb_distance = -log(nb$weight), nearest_weight = nb$weight)
        }
      }
    }, error = function(e) {
      warning("bfb_detect_vectors: skipping vector [", paste(v, collapse = ","),
              "] -- ", conditionMessage(e), call. = FALSE)
      list(admits_bfb = NA, bfb_distance = NA_real_, nearest_weight = NA_real_)
    })
    cache[[key]] <- entry
  }

  na_entry <- list(admits_bfb = NA, bfb_distance = NA_real_, nearest_weight = NA_real_)
  n_unique <- length(unique_keys)

  results <- lapply(seq_along(vec_list), function(i) {
    v <- vec_list[[i]]
    n_seg <- if (is.null(v)) NA_integer_ else length(v)
    max_c <- if (is.null(v)) NA_integer_ else as.integer(max(v))
    key <- keys[i]
    entry <- if (is.na(key)) na_entry else (cache[[key]] %||% na_entry)
    list(n_segments = n_seg, max_count = max_c, admits_bfb = entry$admits_bfb,
         bfb_distance = entry$bfb_distance, nearest_weight = entry$nearest_weight)
  })

  dplyr::tibble(
    n_segments     = vapply(results, `[[`, integer(1), "n_segments"),
    max_count      = vapply(results, `[[`, integer(1), "max_count"),
    admits_bfb     = vapply(results, `[[`, logical(1), "admits_bfb"),
    bfb_distance   = vapply(results, `[[`, numeric(1), "bfb_distance"),
    nearest_weight = vapply(results, `[[`, numeric(1), "nearest_weight"),
    n_unique       = n_unique
  )
}

#' Split per-chromosome copy-number data into p-arm and q-arm
#'
#' Adds an \code{arm} column (\code{"p"}/\code{"q"}) based on each row's
#' position relative to the centromere. This matters for BFB
#' count-vector extraction: genomic coordinates increase from the
#' p-arm's own telomere, through the centromere, to the q-arm's own
#' telomere -- so ascending coordinate order is already
#' telomere-to-centromere for the p-arm, but is centromere-to-telomere
#' (backwards) for the q-arm. Only the q-arm needs reversing to bring
#' both into the telomere-to-centromere convention \code{bfbtools}
#' expects; see \code{.extract_count_vector}'s \code{reverse}
#' argument, which \code{\link{bfb_detect_batch}} sets automatically
#' per arm once you supply a \code{genome}.
#'
#' Centromere positions come from the \code{\link{centromeres}} dataset
#' bundled with this package (real hg19/hg38 coordinates for chromosomes
#' 1-22, X, Y) -- there is no manual-override or fallback-heuristic path,
#' so arm-splitting is always based on real genomic coordinates for a
#' real genome build.
#'
#' @param cna_data a copy-number tibble with (at least) chr/start/end
#'   columns.
#' @param genome genome build to use for centromere positions: one of
#'   \code{"hg19"} (GRCh37) or \code{"hg38"} (GRCh38). Required -- see
#'   \code{\link{centromeres}}.
#' @param chr_col,start_col,end_col column names in \code{cna_data}.
#' @return \code{cna_data} with an added \code{arm} column (\code{"p"}
#'   or \code{"q"}).
#' @export
split_chr_arms <- function(cna_data, genome, chr_col = "chr",
                            start_col = "start", end_col = "end") {
  stopifnot(is.data.frame(cna_data))
  genome <- match.arg(genome, choices = c("hg19", "hg38"))
  centromere <- .centromere_vector_for(genome)

  cna_data$arm <- NA_character_
  chrs <- unique(cna_data[[chr_col]])
  for (chr in chrs) {
    key <- .strip_chr_prefix(as.character(chr))
    if (!key %in% names(centromere)) {
      stop("No centromere position available for chromosome '", key,
           "' under genome = '", genome, "'. bridges bundles centromere data ",
           "for chromosomes 1-22, X, Y only (see ?centromeres).", call. = FALSE)
    }
    idx <- cna_data[[chr_col]] == chr
    cna_data$arm[idx] <- ifelse(cna_data[[start_col]][idx] < centromere[[key]], "p", "q")
  }
  cna_data
}

.strip_chr_prefix <- function(x) sub("^chr", "", x)

#' Extract a BFB count-vector from one row of a copy-number matrix
#'
#' Run-length-encodes a copy-number row, drops zero runs, and treats a
#' single remaining run (no detectable variation) as uninformative.
#'
#' @param reverse whether to reverse the row before RLE-encoding.
#'   Matrix columns are assumed to run in ascending genomic-coordinate
#'   order (as produced by \code{tibble_to_matrix}). Ascending
#'   order is already telomere-to-centromere for a p-arm
#'   (\code{reverse = FALSE}), but is centromere-to-telomere for a
#'   q-arm and needs reversing (\code{reverse = TRUE}) to match
#'   \code{bfbtools}'s expected convention. Defaults to \code{TRUE}
#'   (matching this function's original whole-chromosome behavior) when
#'   arm information isn't available -- see the \code{genome}
#'   argument on \code{\link{bfb_detect_batch}} to enable proper
#'   per-arm handling instead of relying on this default.
#' @noRd
.extract_count_vector <- function(mat, i, reverse = TRUE) {
  row <- mat[i, ]
  if (reverse) row <- rev(row)
  rle_vec <- rle(row)
  v <- rle_vec$values
  v <- v[v != 0]
  if (length(v) <= 1) return(NULL)
  as.integer(v)
}

#' Run BFB detection across a copy-number-data tibble
#'
#' Extracts per-cell, per-allele (and per-chromosome, if \code{cna_data}
#' has more than one) copy-number count-vectors from a
#' \code{bridges}-style \code{cna_data} tibble, and runs
#' \code{\link{bfb_detect_vectors}} on all of them at once (so exact
#' duplicate profiles across cells/alleles are only computed once).
#'
#' ARM AWARENESS: if \code{genome} is supplied (anything other than the
#' default \code{NULL}), each chromosome is first split into p-arm and
#' q-arm via \code{\link{split_chr_arms}} (using real centromere
#' positions for that genome build, from the bundled
#' \code{\link{centromeres}} dataset), and each arm is extracted and
#' detected \emph{separately}, with the correct reversal direction
#' applied automatically per arm (p: no reversal, q: reversed -- see
#' \code{\link{split_chr_arms}}'s documentation for why). The output
#' then has one row per (cell_id, chr, arm, allele) instead of one row
#' per (cell_id, chr, allele), with an added \code{arm} column.
#'
#' If \code{genome} is left \code{NULL} (the default), the entire
#' chromosome is extracted as a single unit with \code{reverse = TRUE}
#' -- this matches this function's original behavior, but note that
#' this is only correct if your data's coordinate convention actually
#' warrants reversing the whole chromosome (e.g. if it only ever
#' represents a single arm to begin with). \strong{If your simulations
#' or data can contain material from both arms of a chromosome, you
#' should supply \code{genome} rather than relying on this default} --
#' see \code{\link{split_chr_arms}}.
#'
#' ASSUMPTIONS (please verify against your actual bridges output and
#' adjust arguments as needed -- I don't have the bridges package
#' available to test this wrapper directly against real data, unlike
#' \code{bfb_detect_vectors()} above, which is fully tested against
#' bfbtools; I have verified the logic below end-to-end, including the
#' arm-splitting/reversal logic, against a minimal stand-in for
#' \code{tibble_to_matrix()} with the assumed interface, but not against
#' the real bridges implementation):
#' \itemize{
#'   \item \code{cna_data} is compatible with
#'     \code{tibble_to_matrix(cna_data, value_column = allele)}
#'     for each allele in \code{alleles}, returning a matrix with one row
#'     per cell and \code{rownames()} giving cell IDs, and with columns
#'     ordered by ascending \code{start}/\code{end} (this is what the
#'     real \code{tibble_to_matrix} does, confirmed from its
#'     source).
#'   \item If \code{cna_data} contains data for more than one
#'     chromosome, there is a column (named by \code{chr_col}, default
#'     \code{"chr"}) identifying it. If no such column is found,
#'     \code{cna_data} is treated as a single chromosome's worth of
#'     data.
#' }
#'
#' @param cna_data a copy-number tibble, e.g. \code{sim$cna_data} from a
#'   \code{bridge_sim()} result.
#' @param alleles character vector of allele labels to process, matching
#'   the \code{value_column} argument to
#'   \code{tibble_to_matrix()} (default \code{c("A", "B")}).
#' @param chr_col name of the chromosome column in \code{cna_data}, if
#'   present (default \code{"chr"}); set to \code{NULL} to force
#'   single-chromosome handling even if such a column exists.
#' @param genome genome build to use for centromere positions, passed to
#'   \code{\link{split_chr_arms}} to enable per-arm handling; one of
#'   \code{"hg19"} or \code{"hg38"}. See "ARM AWARENESS" above. Default
#'   \code{NULL} (no arm-splitting).
#' @param model,min_weight,min_max_count passed to
#'   \code{\link{bfb_detect_vectors}} (see that function for the
#'   trivial-vector filter rationale).
#' @return a tibble with one row per (cell_id, chr, allele) -- or per
#'   (cell_id, chr, arm, allele) if \code{genome} is supplied:
#'   n_segments, max_count, admits_bfb, bfb_distance, nearest_weight,
#'   n_unique (per group).
#' @export
bfb_detect_batch <- function(cna_data, alleles = c("A", "B"), chr_col = "chr",
                              genome = NULL,
                              model = c("poisson", "none"), min_weight = 1e-6,
                              min_max_count = 2) {
  model <- match.arg(model)
  if (!is.data.frame(cna_data)) {
    stop("cna_data must be a data frame/tibble (e.g. sim$cna_data from ",
         "a bridge_sim() result), got: ", class(cna_data)[1], call. = FALSE)
  }

  arm_aware <- !is.null(genome)
  if (arm_aware) {
    cna_data <- split_chr_arms(cna_data, genome = genome, chr_col = chr_col %||% "chr")
  }

  has_chr_col <- !is.null(chr_col) && chr_col %in% names(cna_data)
  chrs <- if (has_chr_col) sort(unique(cna_data[[chr_col]])) else NA_character_
  arms <- if (arm_aware) c("p", "q") else NA_character_

  all_results <- vector("list", length(chrs) * length(arms) * length(alleles))
  idx <- 1L

  for (chr in chrs) {
    cna_chr <- if (has_chr_col) cna_data[cna_data[[chr_col]] == chr, ] else cna_data

    for (arm in arms) {
      cna_sub <- if (arm_aware) cna_chr[cna_chr$arm == arm, ] else cna_chr
      reverse_this <- if (arm_aware) (arm == "q") else TRUE  # see .extract_count_vector docs

      for (allele in alleles) {
        mat <- tibble_to_matrix(cna_sub, value_column = allele)
        cell_ids <- rownames(mat)
        if (is.null(cell_ids)) {
          stop("tibble_to_matrix() returned a matrix with no rownames -- ",
               "can't recover cell_id. See the Assumptions note in ?bfb_detect_batch.",
               call. = FALSE)
        }

        vec_list <- lapply(seq_len(nrow(mat)), .extract_count_vector, mat = mat,
                            reverse = reverse_this)
        res <- bfb_detect_vectors(vec_list, model = model, min_weight = min_weight,
                                   min_max_count = min_max_count)
        res$cell_id <- cell_ids
        res$chr <- if (has_chr_col) chr else NA_character_
        if (arm_aware) res$arm <- arm
        res$allele <- allele

        all_results[[idx]] <- res
        idx <- idx + 1L
      }
    }
  }

  out <- dplyr::bind_rows(all_results)
  id_cols <- if (arm_aware) c("cell_id", "chr", "arm", "allele") else c("cell_id", "chr", "allele")
  dplyr::select(out, dplyr::all_of(id_cols), dplyr::everything())
}
