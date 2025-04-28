
plot_reconstructed_tree = function(x, best_tree, breakpoints, tree_width = 3, raster_quality = 5, use_raster = TRUE) {
  # Retrieve results for best tree

  if (is.list(x)) {
    X = cells2mat(x$cells, x$input_parameters$initial_sequence_length, order = F)
  } else {
    X = x
  }

  L = ncol(X)
  leaves = extract_leaf_vecs_with_names(best_tree)
  lambdas = get_proposed_clusters(leaves, breakpoints, L)

  results <- mixture_of_poissons_cpp(X + 1, lambdas = lambdas + 1)
  leaves_assignment = rownames(lambdas)[results$assignments]

  tree_df = build_phylogeny_df(best_tree)
  tree_df = tree_df %>%
    dplyr::rename(cell_id = seq_name, parent_id = parent_seq) %>%
    dplyr::bind_rows(dplyr::tibble(cell_id = rownames(X), parent_id = leaves_assignment)) %>%
    dplyr::mutate(birth_time = NA, bfb_event = F, hotspot_gained = F)

  tree_data = get_tree_data(list(cell_history = tree_df), branch_lengths = F, full_tree = TRUE)
  tree_ggplot = plot_tree(tree_data, annotate_tip = F)
  ordered_cell_ids <- get_ordered_cell_ids(tree_ggplot$data)

  make_corrupt_tree_heatmap_v2 <- function(tree_ggplot, tree_width, ...) {
    tree_annot_func <- ComplexHeatmap::AnnotationFunction(
      fun = function(index) {
        pushViewport(viewport(height = 1))
        grid.draw(ggplot2::ggplotGrob(tree_ggplot)$grobs[[6]])
        popViewport()
      },
      var_import = list(tree_ggplot = tree_ggplot),
      width = grid::unit(tree_width, "cm"),
      which = "row"
    )
    tree_annot <- ComplexHeatmap::HeatmapAnnotation(
      tree = tree_annot_func, which = "row", show_annotation_name = FALSE
    )

    n_cells <- sum(tree_ggplot$data$isTip)
    tree_hm <- ComplexHeatmap::Heatmap(matrix(ncol = 0, nrow = n_cells), left_annotation = tree_annot, ...)

    return(tree_hm)
  }

  tree_hm = make_corrupt_tree_heatmap_v2(tree_ggplot, tree_width = tree_width)

  clusters = matrix(leaves_assignment, ncol = 1)
  rownames(clusters) = rownames(X)

  mat = X
  high_cn_mask = mat > 10
  mat_row_names = rownames(mat)
  mat = matrix(as.character(mat), ncol = ncol(mat), nrow = nrow(mat))
  mat[high_cn_mask] = "11+"
  mat = data.frame(mat)
  rownames(mat) = mat_row_names
  colvals <- get_colors("CN")

  mat = mat[ordered_cell_ids, , drop = FALSE]
  not_observed_clsuters = ordered_cell_ids[!ordered_cell_ids %in% rownames(clusters)]
  if (length(not_observed_clsuters)) {
    not_observed_clsuters = matrix(NA, nrow = length(not_observed_clsuters), dimnames = list(not_observed_clsuters, NULL))
    clusters = rbind(clusters, not_observed_clsuters)
  }
  clusters = clusters[ordered_cell_ids, , drop=FALSE]

  make_clone_palette <- function(levels) {
    levels <- as.character(levels)  # Ensure levels are characters

    # Assign colors based on the number of levels
    if (length(levels) <= 8) {
      pal <- RColorBrewer::brewer.pal(max(length(levels), 3), "Paired")
    } else {
      pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(max(length(levels), 3), "Paired"))(length(levels))
    }

    # Ensure names are correctly assigned
    names(pal) <- levels
    pal <- pal[levels]

    # Assign a default color for NA values
    if (any(is.na(levels))) {
      pal[is.na(names(pal))] <- "grey"  # Assign grey color for NA clusters
      names(pal)[is.na(names(pal))] <- "NA"
    }

    return(pal)
  }

  # Define annotation colors correctly
  col_anno = make_clone_palette(unique(clusters))

  # Ensure NA clusters are mapped correctly
  if ("NA" %in% names(col_anno)) {
    col_anno["NA"] <- "grey"  # Explicitly set the NA color
  }

  # Create annotation
  cluster_annotation = ComplexHeatmap::rowAnnotation(
    Cluster = clusters,
    col = list(Cluster = col_anno)  # Ensure col is a named list of vectors
  )

  # Heatmap
  copynumber_hm = ComplexHeatmap::Heatmap(
    name = "Copy Number",
    as.matrix(mat),
    col = colvals,
    show_column_names = FALSE,
    show_row_names = FALSE,
    use_raster = use_raster,
    raster_quality = raster_quality,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    left_annotation = cluster_annotation
  )

  # Combine heatmaps
  p = tree_hm + copynumber_hm
  p
}

# Algorithm
# Helper function to find all expandable nodes that contain a given breakpoint
# Helper function: Find all expandable nodes containing a breakpoint
find_expandible_with_bp <- function(node, bp_name, bp_name_rev) {
  if (is.null(node)) return(NULL)

  expandibles <- list()
  seq <- unlist(strsplit(node$seq, ""))

  # Check if current node is expandable and contains breakpoint
  if (is.null(node$left) && is.null(node$right)) {
    for (i in 1:(length(seq) - 1)) {
      pair <- paste0(seq[i:(i + 1)], collapse = "")
      if (pair %in% c(bp_name, bp_name_rev)) {
        expandibles <- c(expandibles, list(list(node = node, cut_idx = i)))
      }
    }
  }

  # Recursively check children (even if current node was expanded)
  c(expandibles,
    find_expandible_with_bp(node$left, bp_name, bp_name_rev),
    find_expandible_with_bp(node$right, bp_name, bp_name_rev))
}

# Modified cut_seq to handle edge cases
cut_seq <- function(seq, cut_idx) {
  fused_seq = c(seq, rev(seq))
  left_seq <- fused_seq[1:cut_idx]
  right_seq <- rev(fused_seq[(cut_idx + 1):length(fused_seq)])

  list(
    left_seq = paste0(left_seq, collapse = ""),
    right_seq = paste0(right_seq, collapse = "")
  )
}

propose_tree <- function(root, bps, N_segments) {
  # If no breakpoints, return the tree unchanged and empty breakpoints
  if (length(bps) == 0) {
    return(list(tree = root, remaining_bps = bps))
  }

  # Select a random breakpoint
  bp_idx <- if (length(bps) > 1) sample(1:length(bps), 1) else 1
  bp <- bps[bp_idx]
  bp_name <- names(bp)
  bp_name_rev <- paste0(rev(strsplit(bp_name, "")[[1]]), collapse = "")

  # Find expandable nodes for this breakpoint
  expandible_options <- find_expandible_with_bp(root, bp_name, bp_name_rev)

  # If no expandable options, return original tree and remove this bp
  if (length(expandible_options) == 0) {
    remaining_bps <- bps[-bp_idx]
    return(list(tree = root, remaining_bps = remaining_bps))
  }

  # Select one random expandable option
  selected <- if (length(expandible_options) > 1) {
    sample(expandible_options, 1)[[1]]
  } else {
    expandible_options[[1]]
  }

  # Apply the breakpoint to create a new tree
  new_root <- apply_breakpoint(root, selected$node, selected$cut_idx, N_segments)

  # Remove the used breakpoint from the list
  remaining_bps <- bps[-bp_idx]

  # Return both the new tree and remaining breakpoints
  return(list(tree = new_root, remaining_bps = remaining_bps))
}

# Main tree building with proper updates
build_tree <- function(root, bps) {
  while (length(bps) > 0) {
    # Select breakpoint
    bp <- if (length(bps) > 1) sample(bps, 1) else bps[1]
    bp_name <- names(bp)
    bp_name_rev <- paste0(rev(strsplit(bp_name, "")[[1]]), collapse = "")

    # Find expandable nodes
    expandible_options <- find_expandible_with_bp(root, bp_name, bp_name_rev)

    if (length(expandible_options) == 0) {
      bps <- bps[names(bps) != bp_name]
      next
    }

    # Select random option
    selected <- if (length(expandible_options) > 1) {
      sample(expandible_options, 1)[[1]]
    } else {
      expandible_options[[1]]
    }

    # Apply breakpoint (returns NEW modified tree)
    root <- apply_breakpoint(root, selected$node, selected$cut_idx, N)
    bps <- bps[names(bps) != bp_name]
  }
  return(root)
}

# Applies breakpoint and returns updated tree
apply_breakpoint <- function(current_node, target_node, cut_idx, N) {
  if (identical(current_node, target_node)) {
    seq <- unlist(strsplit(current_node$seq, ""))
    seqs <- cut_seq(seq, cut_idx)
    current_node$left <- init_cell(seqs$left_seq, N)
    current_node$right <- init_cell(seqs$right_seq, N)
    return(current_node)
  }

  if (!is.null(current_node$left)) {
    current_node$left <- apply_breakpoint(current_node$left, target_node, cut_idx, N)
  }
  if (!is.null(current_node$right)) {
    current_node$right <- apply_breakpoint(current_node$right, target_node, cut_idx, N)
  }
  return(current_node)
}

# Install if needed
library(data.tree)

# Convert your tree to a data.tree structure
convert_to_dataTree <- function(node, parent = NULL, convert_name = FALSE) {
  if (is.null(node)) return(NULL)

  if (convert_name) {
    name = paste0("[",paste0(node$vec, collapse = ","), "]")
  } else {
    name = node$seq
  }
  curr <- Node$new(name)
  if (!is.null(parent)) parent$AddChildNode(curr)

  convert_to_dataTree(node$left, curr, convert_name)
  convert_to_dataTree(node$right, curr, convert_name)

  return(curr)
}

extract_leaf_vecs_with_names <- function(node) {
  if (is.null(node)) return(NULL)

  if (is.null(node$left) && is.null(node$right)) {
    return(stats::setNames(list(node$vec), node$seq))
  }

  c(extract_leaf_vecs_with_names(node$left),
    extract_leaf_vecs_with_names(node$right))
}

seq2counts = function(seq, N) {
  observed = table(unlist(strsplit(seq, "")))
  vec = rep(0, N)
  names(vec) = LETTERS[1:N]
  for (n in names(observed)) {
    vec[n] = observed[n]
  }
  vec
}

init_cell = function(seq, N) {
  seq_vec = seq2counts(seq, N)
  list(seq = seq, vec = seq_vec, left = NULL, right = NULL)
}


get_proposed_clusters = function(leaves, breakpoints, L) {
  proposed_clusters = lapply(leaves, function(l) {
    v = rep(0, L)
    expanded_bp = c(0, breakpoints, L)
    for (idx in 1:(length(breakpoints)+1)) {
      limits = expanded_bp[idx:(idx+1)]
      v[(limits[1]+1):limits[2]] = l[idx]
    }
    v
  }) %>% do.call("rbind", .)
  proposed_clusters
}


## Mixture of Poissons - Parameter Estimation with Fixed Lambda Values
# Input:
# - X: N x L matrix where each row is an observation
# - lambdas: K x L matrix where each row is a lambda vector for one mixture component
# - max_iter: Maximum number of EM iterations
# - tol: Convergence tolerance

mixture_of_poissons <- function(X, lambdas, max_iter = 100, tol = 1e-6) {
  N <- nrow(X)
  L <- ncol(X)
  K <- nrow(lambdas)

  # Initialize mixture weights uniformly
  pi_k <- rep(1/K, K)

  # Initialize log-likelihood
  prev_ll <- -Inf

  # Store cluster assignments
  assignments <- rep(NA, N)

  # EM algorithm
  for (iter in 1:max_iter) {
    # E-step: Calculate responsibilities (posterior probabilities)
    responsibilities <- matrix(0, N, K)

    for (i in 1:N) {
      for (k in 1:K) {
        # Calculate Poisson probabilities for observation i under component k
        log_probs <- stats::dpois(X[i,], lambda = lambdas[k,], log = TRUE)
        # Sum log probabilities (equivalent to product of probabilities)
        log_prob_i_k <- sum(log_probs)
        # Store log(pi_k * P(x_i|lambda_k))
        responsibilities[i, k] <- log(pi_k[k]) + log_prob_i_k
      }

      # Normalize to get probabilities (with numerical stability)
      max_resp <- max(responsibilities[i,])
      responsibilities[i,] <- exp(responsibilities[i,] - max_resp)
      responsibilities[i,] <- responsibilities[i,] / sum(responsibilities[i,])
    }

    # M-step: Update mixture weights
    pi_k <- colMeans(responsibilities)

    # Calculate log-likelihood
    ll <- 0
    for (i in 1:N) {
      log_probs_i <- numeric(K)
      for (k in 1:K) {
        component_probs <- stats::dpois(X[i,], lambda = lambdas[k,], log = TRUE)
        log_probs_i[k] <- log(pi_k[k]) + sum(component_probs)
      }

      # Log-sum-exp trick for numerical stability
      max_log_prob <- max(log_probs_i)
      ll <- ll + max_log_prob + log(sum(exp(log_probs_i - max_log_prob)))
    }

    # Check convergence
    if (abs(ll - prev_ll) < tol) {
      break
    }

    prev_ll <- ll
  }

  # Determine cluster assignments (hard assignment)
  assignments <- apply(responsibilities, 1, which.max)

  # Return results
  return(list(
    pi_k = pi_k,                      # Estimated mixture weights
    responsibilities = responsibilities, # Soft assignments (probabilities)
    assignments = assignments,         # Hard assignments
    log_likelihood = prev_ll           # Final log-likelihood
  ))
}

# Function to traverse the tree and extract sequence relationships
build_phylogeny_df <- function(node, parent_seq = NA) {
  if (is.null(node)) {
    return(NULL)
  }

  # Extract current node's sequence
  current_seq <- node$seq

  # Create data frame row for current node
  current_df <- data.frame(
    seq_name = current_seq,
    parent_seq = parent_seq,
    stringsAsFactors = FALSE
  )

  # Recursively process left and right children
  left_df <- build_phylogeny_df(node$left, current_seq)
  right_df <- build_phylogeny_df(node$right, current_seq)

  # Combine results
  result_df <- rbind(current_df, left_df, right_df)
  return(result_df)
}

# Function to recursively build Newick string from tree structure
build_newick_string <- function(node) {
  if (is.null(node)) {
    return("")
  }

  # Get current node's sequence
  current_seq <- node$seq

  # If this is a leaf node (no children)
  if (is.null(node$left) && is.null(node$right)) {
    return(current_seq)
  }

  # Process children
  left_newick <- build_newick_string(node$left)
  right_newick <- build_newick_string(node$right)

  # Build the current node's Newick string
  children_newick <- ""

  # Add left child if it exists
  if (left_newick != "") {
    children_newick <- left_newick
  }

  # Add right child if it exists, with comma if needed
  if (right_newick != "") {
    if (children_newick != "") {
      children_newick <- paste0(children_newick, ",", right_newick)
    } else {
      children_newick <- right_newick
    }
  }

  # Construct final Newick string for this subtree
  if (children_newick != "") {
    return(paste0("(", children_newick, ")", current_seq))
  } else {
    return(current_seq)
  }
}

# Get the full Newick string and add the terminating semicolon
get_newick_tree <- function(tree) {
  newick <- build_newick_string(tree)
  return(paste0(newick, ";"))
}
