#' Test for BFB signal along the branches of an inferred phylogeny
#'
#' For every chromosome and allele in a \code{\link{fit}} result, counts the
#' internal nodes of the tree whose children are better explained by a BFB
#' event than by greedy gain/loss steps, and tests whether that count is
#' larger than expected under a background rate \code{threshold}.
#'
#' This test counts BFB-like nodes over the whole tree, so BFB confined to a
#' part of the tree can be diluted; \code{\link{detect_bfb_clades}} scans every
#' clade instead.
#'
#' This is one of two complementary BFB-detection approaches in the package.
#' \code{detect_bfb_branches()} works on a fitted tree; \code{\link{bfb_detect_batch}}
#' works directly on each cell's copy-number profile, with no tree needed,
#' and asks whether that profile is compatible with any BFB schedule.
#'
#' @param fit A result of \code{\link{fit}}.
#' @param threshold Background probability that a branch looks BFB-like by
#'   chance; the null proportion of the one-sided binomial test.
#'   Default: 0.005
#'
#' @return A data frame with one row per chromosome-allele combination and
#'   columns:
#' \describe{
#'   \item{mean}{Fraction of internal nodes flagged as BFB-like.}
#'   \item{N_trials}{Number of internal nodes tested.}
#'   \item{N_successes}{Number of internal nodes flagged as BFB-like.}
#'   \item{chr}{Chromosome identifier.}
#'   \item{allele}{Allele.}
#'   \item{p.value}{One-sided binomial test p-value.}
#'   \item{adj.pval}{Benjamini-Hochberg adjusted p-value across all
#'     chromosome-allele combinations.}
#'   \item{threshold}{The \code{threshold} used.}
#' }
#'
#' @details
#' During \code{fit()}, ancestral profiles are reconstructed bottom-up along
#' the tree. At each internal node, the two children are compared with a
#' greedy gain/loss cost G and a BFB cost B; \eqn{\delta = G - B}. A node with
#' \eqn{\delta > 0} is counted as BFB-like. For each chromosome-allele, the
#' number of BFB-like nodes is tested against
#' \code{binom.test(k, N, p = threshold, alternative = "greater")}, and
#' p-values are adjusted with \code{stats::p.adjust(method = "BH")}.
#'
#' @examples
#' \dontrun{
#' sim <- bridge_sim(chromosomes = "8", bfb_allele = "8:A",
#'                   max_cells = 128, lambda = 2)
#' res <- fit(data = sim$cna_data, alleles = c("A", "B"))
#' bfb_df <- detect_bfb_branches(res, threshold = 0.005)
#' head(bfb_df)
#' # Chromosome 8 allele A should show high mean and low adj.pval
#' }
#' @export
detect_bfb_branches = function(fit, threshold = .005) {
  alleles = names(fit$all_input_Xs[[1]])
  chromosomes = names(fit$all_input_Xs)

  pseudo_cells_test_df = lapply(alleles, function(allele) {
    lapply(chromosomes, function(chr) {
      pseudo_cell_bin_test(
        fit = fit,
        chr = chr,
        allele = allele,
        threshold = threshold
      )
    }) %>%
      do.call(dplyr::bind_rows, .)
  }) %>%
    do.call(dplyr::bind_rows, .)

  pseudo_cells_test_df$adj.pval = stats::p.adjust(pseudo_cells_test_df$p.value, method = "BH")
  pseudo_cells_test_df
}

#' Perform Binomial Test on Pseudo-Cell Delta Values
#'
#' This function reconstructs a phylogenetic tree and performs a binomial test
#' to assess whether the proportion of positive delta values exceeds a given threshold.
#' Delta values represent the difference between G distance and B distance costs.
#'
#' @param fit A fitted model object containing tree structure and distance functions
#' @param chr Character string specifying the chromosome
#' @param allele Character string specifying the allele type
#' @param threshold Numeric threshold probability for the binomial test
#'
#' @return A tibble containing:
#'   \itemize{
#'     \item mean: Estimated proportion of positive delta values
#'     \item N_trials: Total number of delta values
#'     \item N_successes: Number of positive delta values
#'     \item chr: Chromosome identifier
#'     \item allele: Allele type
#'     \item p.value: P-value from binomial test
#'     \item threshold: Threshold used for testing
#'   }
#'
#' @details The binomial test uses alternative = "greater" to test if the proportion
#' of positive deltas significantly exceeds the threshold.
#' @keywords internal
pseudo_cell_bin_test = function(fit, chr, allele, threshold) {
  r = fit$reconstructions[[chr]][[allele]]  # use precomputed
  if (is.null(r)) r = reconstruct_tree(fit, chr, allele)
  delta_values <- unlist(r$deltas)
  delta_values <- delta_values[!is.na(delta_values)]  # NA: WGD nodes, not tested
  N = length(delta_values)
  k = sum(delta_values > 0)

  test = stats::binom.test(k, N, p = threshold, alternative = "greater")
  dplyr::tibble(
    mean = test$estimate,
    N_trials = N,
    N_successes = k,
    chr = chr,
    allele = allele,
    p.value = test$p.value,
    threshold = threshold
  )
}

get_post_order_internal_nodes = function(tree) {
  Ntip <- length(tree$tip.label)
  Nnode <- tree$Nnode

  # Build children list
  children_list <- vector("list", max(tree$edge))
  for (i in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    children_list[[parent]] <- c(children_list[[parent]], child)
  }

  # Iterative post-order traversal using a stack
  stack <- list(list(node = Ntip + 1, visited = FALSE))
  post_order <- c()

  while (length(stack) > 0) {
    current <- stack[[length(stack)]]
    stack <- stack[-length(stack)]  # pop

    if (current$visited) {
      if (current$node > Ntip) {
        post_order <- c(post_order, current$node)
      }
    } else {
      # Mark as visited and push back
      stack <- c(stack, list(list(node = current$node, visited = TRUE)))

      # Push children (in reverse order for correct post-order)
      children <- children_list[[current$node]]
      if (!is.null(children)) {
        internal_children <- children[children > Ntip]
        for (child in rev(internal_children)) {
          stack <- c(stack, list(list(node = child, visited = FALSE)))
        }
      }
    }
  }
  post_order
}

compute_reconstructions = function(fit, chromosomes = NULL, alleles = NULL) {
  if (is.null(chromosomes)) {
    chromosomes <- names(fit$all_input_Xs)
  } else {
    chromosomes <- intersect(chromosomes, names(fit$all_input_Xs))
    if (length(chromosomes) == 0) {
      stop("None of the specified chromosomes are present in the fit object.")
    }
  }

  if (is.null(alleles)) {
    alleles <- names(fit$all_input_Xs[[chromosomes[1]]])
  } else {
    # Ensure all selected chromosomes have the requested alleles
    common_alleles <- Reduce(intersect, lapply(fit$all_input_Xs[chromosomes], names))
    alleles <- intersect(alleles, common_alleles)
    if (length(alleles) == 0) {
      stop("None of the specified alleles are present in the selected chromosomes of the fit object.")
    }
  }

  new_edges = NULL
  whole_history = dplyr::tibble()
  reconstructions = list()

  for (chr in chromosomes) {
    reconstructions[[chr]] <- list()
    for (all in alleles) {
      reconstruction = reconstruct_tree(fit = fit, chr = chr, allele = all)
      reconstructions[[chr]][[all]] <- reconstruction
      whole_history = dplyr::bind_rows(whole_history, reconstruction$history)
      if (!is.null(new_edges)) {
        new_edges = new_edges + reconstruction$new_edge_lengths
      } else {
        new_edges = reconstruction$new_edge_lengths
      }
    }
  }

  # Each WGD edge costs wgd_penalty per WGD, once for the whole genome
  if (!is.null(fit$wgd)) {
    dk <- fit$wgd$k_node[fit$tree$edge[, 2]] - fit$wgd$k_node[fit$tree$edge[, 1]]
    new_edges <- new_edges + fit$wgd$penalty * dk
  }

  list(
    reconstructions = reconstructions,
    edge.length = new_edges,
    history = whole_history
  )
}

#' Reconstruct Phylogenetic Tree with Copy Number Analysis
#'
#' This function reconstructs a phylogenetic tree by traversing internal nodes and
#' computing distances between cell profiles. For each internal node, it calculates
#' both B distance (breakpoint-aware) and G distance (general) between child nodes,
#' then selects the reconstruction method that minimizes cost.
#'
#' @param fit A fitted model object containing:
#'   \itemize{
#'     \item all_input_Xs: List of input data matrices by chromosome and allele
#'     \item b_dist_func: B distance function for breakpoint-aware reconstruction
#'     \item g_dist_func: G distance function for general reconstruction
#'     \item tree: Phylogenetic tree structure (ape format)
#'   }
#' @param chr Character string specifying the chromosome to analyze
#' @param allele Character string specifying the allele type ("CN" for copy number, others for allelic)
#'
#' @return A list containing:
#'   \itemize{
#'     \item internal_nodes: List of reconstructed internal node profiles
#'     \item deltas: List of delta values (G cost - B cost) for each internal node
#'     \item merged_profiles: List of merged profiles when B distance is chosen
#'     \item final_input: Final working input matrix after all reconstructions
#'     \item processed_nodes: Number of internal nodes processed
#'     \item chr: Chromosome identifier
#'     \item allele: Allele type
#'   }
#'
#' @details The function:
#' \enumerate{
#'   \item Sets target values based on allele type (2 for "CN", 1 for others)
#'   \item Processes internal nodes in post-order traversal
#'   \item For each internal node, computes B and G distances between children
#'   \item Calculates delta = G_cost - B_cost
#'   \item Chooses B reconstruction if delta > 0, otherwise G reconstruction
#'   \item Creates pseudo-cells for internal nodes and updates working dataset
#' }
#'
#' @note The tree edge lengths are set to 1 for uniform weighting. Node labels
#' are automatically generated if not present.
#' @param store_profiles Logical; whether to retain each internal node's
#'   reconstructed profile in the returned \code{internal_nodes} list
#'   (set FALSE to save memory on large trees when only the final
#'   working input/deltas are needed). Default TRUE.
#'
#' @keywords internal
reconstruct_tree <- function(fit, chr, allele, store_profiles = TRUE) {
  # Get distances
  B_dist <- if (!is.null(fit$b_dist_func)) get_b_dist(fit$b_dist_func) else function(a, b, penalty) list(cost = Inf)
  G_dist <- get_g_dist(fit$g_dist_func)

  # Inputs and tree prep
  target_val <- if (allele == "CN") 2 else 1
  input <- fit$all_input_Xs[[chr]][[allele]]  # rows are cells (tips), cols are bins

  tree <- fit$tree
  # WGD counts per node (fit(wgd = TRUE)); numbering refers to fit$tree
  k_node <- if (!is.null(fit$wgd) && ape::is.rooted(tree)) fit$wgd$k_node else NULL
  if (!ape::is.rooted(tree)) tree <- phangorn::midpoint(tree)
  tree$edge.length[] <- 1

  if (is.null(tree$node.label)) tree$node.label <- paste0("N", seq_len(tree$Nnode))

  n_tip  <- length(tree$tip.label)
  n_int  <- tree$Nnode
  n_node <- n_tip + n_int
  p      <- ncol(input)

  # Names by phylo indexing (tips 1..n_tip, internals n_tip+1..n_node)
  all_names <- c(tree$tip.label, tree$node.label)

  # Postorder internal traversal
  post_order_nodes <- get_post_order_internal_nodes(tree) # numeric IDs of internal nodes in post-order

  # Children map (one pass)
  edge <- tree$edge
  children_of <- split(edge[, 2], edge[, 1])  # list: parentID -> integer vector of children (length 2)

  # Preallocate profiles and means
  profiles <- matrix(NA_real_, nrow = n_node, ncol = p,
                     dimnames = list(c(tree$tip.label, paste0("pseudo_", tree$node.label)), colnames(input)))

  # Reorder input to match phylo tip order
  profiles[seq_len(n_tip), ] <- input[tree$tip.label, , drop = FALSE]
  row_means <- rep(NA_real_, n_node)
  row_means[seq_len(n_tip)] <- rowMeans(profiles[seq_len(n_tip), , drop = FALSE])

  # branch costs indexed by child node id
  branch_costs_vec <- rep(NA_real_, n_node)

  # internal node results keyed by node id (as character) for compatibility
  idx_names <- as.character(post_order_nodes)
  internal_nodes <- vector("list", length(post_order_nodes)); names(internal_nodes) <- idx_names
  merged_profiles <- vector("list", length(post_order_nodes)); names(merged_profiles) <- idx_names
  deltas_vec <- stats::setNames(numeric(length(post_order_nodes)), idx_names)

  # history: two rows per internal node
  H <- 2L * length(post_order_nodes)
  hist_cell_id   <- character(H)
  hist_parent_id <- character(H)
  hist_event     <- character(H)
  h <- 0L

  event_type <- function(d) {
    pos <- any(d > 0)
    neg <- any(d < 0)
    if (!pos && !neg) "none" else if (pos && neg) "gain and loss" else if (pos) "gain" else "loss"
  }

  for (node_id in post_order_nodes) {
    kids <- children_of[[as.character(node_id)]]
    # safety: ensure binary
    if (length(kids) != 2L) next
    kL <- kids[1]; kR <- kids[2]

    left  <- profiles[kL, ]
    right <- profiles[kR, ]

    # Compute row means lazily for ordering in B distance
    mL <- row_means[kL]; if (is.na(mL)) { mL <- mean(left);  row_means[kL] <- mL }
    mR <- row_means[kR]; if (is.na(mR)) { mR <- mean(right); row_means[kR] <- mR }

    # WGD counts of this node and its children (all 0 without fit$wgd)
    k_here <- if (is.null(k_node)) 0 else k_node[node_id]
    fL <- if (is.null(k_node)) 1 else 2^(k_node[kL] - k_here)
    fR <- if (is.null(k_node)) 1 else 2^(k_node[kR] - k_here)

    if (fL > 1 || fR > 1) {
      # One or both children descend from a WGD after this node: the ancestor
      # is the best pre-doubling profile, and BFB is not tested here.
      if (fL > 1 && fR > 1) {
        g_h <- G_dist(floor(left / fL), floor(right / fR), target_val = target_val * 2^k_here)
        internal_cell <- g_h$ancestor
        left_cost  <- g_h$a_cost + G_dist(fL * internal_cell, left, target_val = target_val)$cost
        right_cost <- g_h$b_cost + G_dist(fR * internal_cell, right, target_val = target_val)$cost
      } else if (fR > 1) {
        r <- wgd_ancestor_cpp(left, right, fR)
        internal_cell <- r$x
        left_cost  <- r$steps
        right_cost <- r$after
      } else {
        r <- wgd_ancestor_cpp(right, left, fL)
        internal_cell <- r$x
        right_cost <- r$steps
        left_cost  <- r$after
      }
      left_event  <- if (fL > 1) "wgd" else event_type(left - internal_cell)
      right_event <- if (fR > 1) "wgd" else event_type(right - internal_cell)
      left_vec  <- left
      right_vec <- right
      deltas_vec[as.character(node_id)] <- NA_real_
    } else {
      # Greedy distance; inside clades with k WGDs the ancestor is chosen
      # closest to the doubled normal level (target_val * 2^k)
      g_d <- G_dist(left, right, target_val = target_val * 2^k_here)

      # BFB distance (order matters)
      if (mL > mR) {
        b_d <- B_dist(left, right, penalty = 0)
        b_left_cost  <- b_d$a_cost
        b_right_cost <- b_d$b_cost
        b_left_vec   <- b_d$left
        b_right_vec  <- b_d$right
      } else {
        b_d <- B_dist(right, left, penalty = 0)
        b_left_cost  <- b_d$b_cost
        b_right_cost <- b_d$a_cost
        b_left_vec   <- b_d$right
        b_right_vec  <- b_d$left
      }

      delta <- g_d$cost - b_d$cost
      deltas_vec[as.character(node_id)] <- delta

      if (is.finite(b_d$cost) && delta > 0) {
        # BFB wins
        left_cost  <- b_left_cost
        right_cost <- b_right_cost
        internal_cell <- b_d$ancestor
        left_event  <- "bfb"
        right_event <- "bfb"
        left_vec  <- b_left_vec
        right_vec <- b_right_vec
      } else {
        # Greedy wins (or BFB not available)
        left_cost  <- g_d$a_cost
        right_cost <- g_d$b_cost
        internal_cell <- g_d$ancestor
        left_event  <- event_type(g_d$a_delta)
        right_event <- event_type(g_d$b_delta)
        left_vec <- g_d$ancestor + g_d$a_delta
        right_vec <- g_d$ancestor + g_d$b_delta
      }
    }

    # Save branch costs by child id
    branch_costs_vec[kL] <- left_cost
    branch_costs_vec[kR] <- right_cost

    # Save internal node profile
    profiles[node_id, ] <- if (length(internal_cell) == 1L) rep(internal_cell, p) else internal_cell
    row_means[node_id]  <- mean(profiles[node_id, ])

    # Store results (keyed by node id as character)
    internal_nodes[[as.character(node_id)]] <- profiles[node_id, ]
    if (store_profiles) {
      merged_profiles[[as.character(node_id)]] <- list(left = left_vec, right = right_vec)
    } else {
      merged_profiles[[as.character(node_id)]] <- NULL
    }

    # History rows (use same names as your original: leaves keep names; internals are "pseudo_<name>")
    parent_name <- paste0("pseudo_", all_names[node_id])
    childL_name <- if (kL <= n_tip) all_names[kL] else paste0("pseudo_", all_names[kL])
    childR_name <- if (kR <= n_tip) all_names[kR] else paste0("pseudo_", all_names[kR])

    h <- h + 1L
    hist_cell_id[h]   <- childL_name
    hist_parent_id[h] <- parent_name
    hist_event[h]     <- left_event

    h <- h + 1L
    hist_cell_id[h]   <- childR_name
    hist_parent_id[h] <- parent_name
    hist_event[h]     <- right_event
  }

  # Trim preallocated history (in case of non-binary nodes)
  if (h < H) {
    keep <- seq_len(h)
    hist_cell_id   <- hist_cell_id[keep]
    hist_parent_id <- hist_parent_id[keep]
    hist_event     <- hist_event[keep]
  }
  history <- data.frame(
    cell_id = hist_cell_id,
    parent_id = hist_parent_id,
    event = hist_event,
    chr = chr,
    allele = allele,
    stringsAsFactors = FALSE
  )

  # New edge lengths from child-specific costs, default 1 if missing
  new_edge_lengths <- rep(1, nrow(edge))
  bc_child <- branch_costs_vec[edge[, 2]]
  new_edge_lengths[is.finite(bc_child)] <- bc_child[is.finite(bc_child)]

  # Root node & final_input (root profile as 1xP matrix)
  root_id <- setdiff(edge[, 1], edge[, 2])[1]
  final_input <- matrix(profiles[root_id, ], nrow = 1, dimnames = list(paste0("pseudo_", all_names[root_id]), colnames(input)))

  # Convert vectors to lists keyed like your original
  branch_costs <- as.list(branch_costs_vec[is.finite(branch_costs_vec)])
  names(branch_costs) <- as.character(which(is.finite(branch_costs_vec)))

  deltas <- as.list(deltas_vec)

  list(
    internal_nodes = internal_nodes,
    deltas = deltas,
    merged_profiles = merged_profiles,
    final_input = final_input,
    processed_nodes = length(internal_nodes),
    chr = chr,
    allele = allele,
    branch_costs = branch_costs,
    new_edge_lengths = new_edge_lengths,
    history = history
  )
}
