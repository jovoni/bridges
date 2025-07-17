#' Detect Breakpoint-Free Branches (BFB) in Phylogenetic Trees
#'
#' This function performs binomial tests across all chromosomes and alleles to detect
#' breakpoint-free branches in a phylogenetic tree reconstruction. It uses pseudo-cell
#' analysis to identify branches where copy number changes occur without breakpoints.
#'
#' @param fit A fitted model object containing:
#'   \itemize{
#'     \item all_input_Xs: List of chromosome data with allele information
#'     \item b_dist_func: B distance function
#'     \item g_dist_func: G distance function
#'     \item tree: Phylogenetic tree structure
#'   }
#' @param threshold Numeric threshold for binomial test probability (default: 0.005)
#'
#' @return A data frame containing binomial test results for each chromosome-allele combination:
#'   \itemize{
#'     \item mean: Estimated proportion of successes
#'     \item N_trials: Total number of trials
#'     \item N_successes: Number of successes
#'     \item chr: Chromosome identifier
#'     \item allele: Allele type
#'     \item p.value: P-value from binomial test
#'     \item threshold: Threshold used for testing
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'fitted_model' is a properly fitted model object
#' bfb_results <- detect_bfb(fitted_model, threshold = 0.01)
#' }
detect_bfb = function(fit, threshold = .005) {
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

  pseudo_cells_test_df
}

#' Perform T-Test on Pseudo-Cell Delta Values
#'
#' This function reconstructs a phylogenetic tree for a specific chromosome and allele,
#' then performs a one-sample t-test on the delta values (difference between G and B distances).
#'
#' @param fit A fitted model object containing tree structure and distance functions
#' @param chr Character string specifying the chromosome
#' @param allele Character string specifying the allele type
#' @param mu Numeric value for the null hypothesis mean (default: 0)
#' @param sign Logical indicating whether to use sign of delta values (default: FALSE)
#'
#' @return A tibble containing:
#'   \itemize{
#'     \item mean: Sample mean of delta values
#'     \item chr: Chromosome identifier
#'     \item allele: Allele type
#'     \item p.value: P-value from t-test (NA if all values are identical)
#'     \item greater: Logical indicating if mean > mu
#'   }
#'
#' @details If all delta values are identical, no statistical test is performed and
#' p.value is set to NA.
pseudo_cell_t_test = function(
    fit,
    chr,
    allele,
    mu = 0,
    sign = F
) {
  r <- reconstruct_tree(fit, chr, allele)

  delta_values <- unlist(r$deltas)
  if (sign) delta_values <- sign(delta_values)

  if (length(unique(delta_values)) == 1) {
    estimate = mean(delta_values)
    dplyr::tibble(
      mean = estimate,
      chr = chr,
      allele = allele,
      p.value = NA,
      greater = estimate > mu
    )
  } else {
    test = stats::t.test(delta_values, alternative = "greater", mu = mu)
    dplyr::tibble(
      mean = test$estimate,
      chr = chr,
      allele = allele,
      p.value = test$p.value,
      greater = test$estimate > mu
    )
  }
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
pseudo_cell_bin_test = function(fit, chr, allele, threshold) {
  r = fit$reconstructions[[chr]][[allele]]  # use precomputed
  delta_values <- unlist(r$deltas)
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
reconstruct_tree = function(fit, chr, allele) {
  B_dist = get_b_dist(fit$b_dist_func)
  G_dist = get_g_dist(fit$g_dist_func)

  # Set target value based on allele type
  target_val <- ifelse(allele == "CN", 2, 1)

  # Get input data for the specified chromosome and allele
  input <- fit$all_input_Xs[[chr]][[allele]]

  # Get tree structure
  tree <- fit$tree

  if (!ape::is.rooted(tree)) {
    tree = phangorn::midpoint(tree)
  }
  tree$edge.length = rep(1, length(tree$edge.length))

  # Get post-order traversal of internal nodes to ensure children are processed before parents
  post_order_nodes <- get_post_order_internal_nodes(tree)
  if (is.null(tree$node.label)) tree$node.label <- paste0("N", seq_len(tree$Nnode))

  # Initialize results storage
  internal_nodes <- list()
  deltas <- list()
  merged_profiles = list()

  # Initialize branch costs storage
  branch_costs <- list()

  # Create a working copy of input data that will be modified
  working_input <- input

  # Create a mapping from node to its representative cell name
  node_to_cell <- list()

  leaves_names = tree$tip.label
  node_names = tree$node.label
  all_names = c(leaves_names, node_names)

  # Initialize leaf nodes to map to themselves
  for (tip_name in tree$tip.label) {
    node_to_cell[[tip_name]] <- tip_name
  }

  # Process each internal node in post-order
  #plot(tree, show.node.label = T)
  #node_id = 13

  node_id = post_order_nodes[1]
  history = dplyr::tibble()
  for (node_id in post_order_nodes) {
    #print(which(post_order_nodes == node_id))
    # print(working_input)
    # print(all_names[node_id])
    # Get children of this node
    children <- get_node_children(tree, node_id)

    if (length(children) == 2) {
      # Get the cell representatives for left and right children
      left_child <- children[1]
      right_child <- children[2]

      # Get cell names (either original tip names or pseudo-cell names)
      if (left_child <= length(tree$tip.label)) {
        left_cell_name <- all_names[left_child]
      } else {
        left_cell_name <- paste0("pseudo_", all_names[left_child])
      }

      if (right_child <= length(tree$tip.label)) {
        right_cell_name <- all_names[right_child]
      } else {
        right_cell_name <- paste0("pseudo_", all_names[right_child])
      }

      #if ("cell_14" %in% c(left_cell_name, right_cell_name)) stop()

      # Extract left and right cells from working input
      left_cell <- working_input[left_cell_name, ]
      right_cell <- working_input[right_cell_name, ]

      # Calculate G distance
      g_d <- G_dist(left_cell, right_cell, target_val = target_val)

      # Calculate B distance (order matters based on mean)
      if (mean(left_cell) > mean(right_cell)) {
        b_d <- B_dist(left_cell, right_cell, penalty = 0)
      } else {
        b_d <- B_dist(right_cell, left_cell, penalty = 0)
      }

      # Calculate delta
      delta <- g_d$cost - b_d$cost

      # Choose internal node based on which distance is lower
      if (delta > 0) {
        if (mean(left_cell) > mean(right_cell)) {
          left_cost = b_d$a_cost
          right_cost = b_d$b_cost

          merged_profiles[[as.character(node_id)]] = list(
            left = b_d$left,
            right = b_d$right
          )

        } else {
          left_cost = b_d$b_cost
          right_cost = b_d$a_cost

          merged_profiles[[as.character(node_id)]] = list(
            left = b_d$right,
            right = b_d$left
          )
        }

        left_event = right_event = "bfb"

        internal_node_cell <- b_d$ancestor

      } else {
        left_cost = g_d$a_cost
        right_cost = g_d$b_cost

        get_event_type = function(delta) {
          signs = unique(sign(delta))
          if (all(signs == 0)) return("none")
          if ((1 %in% signs) & (-1 %in% signs)) return("gain and loss")
          if ((1 %in% signs)) return("gain")
          if ((-1 %in% signs)) return("loss")
        }

        left_event = get_event_type(g_d$a_delta)
        right_event = get_event_type(g_d$b_delta)

        internal_node_cell <- g_d$ancestor
        merged_profiles[[as.character(node_id)]] = list(
          left = g_d$left,
          right = g_d$right
        )
      }

      # Store branch costs for left and right children
      branch_costs[[as.character(left_child)]] <- left_cost
      branch_costs[[as.character(right_child)]] <- right_cost

      if (length(internal_node_cell) == 1) {
        internal_node_cell = matrix(internal_node_cell, nrow = 1, ncol = 1)
      }

      # Store results
      internal_nodes[[as.character(node_id)]] <- internal_node_cell
      deltas[[as.character(node_id)]] <- delta

      # Add pseudo-cell to working input
      pseudo_cell_name <- paste0("pseudo_", all_names[node_id])
      working_input <- rbind(working_input, internal_node_cell)
      rownames(working_input)[nrow(working_input)] <- pseudo_cell_name

      # Remove the cells that were combined (but keep other cells)
      working_input <- working_input[
        !rownames(working_input) %in% c(left_cell_name, right_cell_name),
        ,
        drop = F
      ]

      history = dplyr::bind_rows(
        history,
        dplyr::tibble(cell_id = left_cell_name, parent_id = pseudo_cell_name, event = left_event),
        dplyr::tibble(cell_id = right_cell_name, parent_id = pseudo_cell_name, event = right_event)
      )

      #print(rownames(working_input))
    }
  }
  history = history %>% dplyr::mutate(chr = chr, allele = allele)

  # Create new branch lengths based on costs
  new_edge_lengths <- rep(NA, nrow(tree$edge))
  for (i in 1:nrow(tree$edge)) {
    child_node <- tree$edge[i, 2]
    if (as.character(child_node) %in% names(branch_costs)) {
      new_edge_lengths[i] <- branch_costs[[as.character(child_node)]]
    } else {
      # If no cost available, keep original length or set to 1
      new_edge_lengths[i] <- 1
    }
  }

  # Return results as a list
  list(
    internal_nodes = internal_nodes,
    deltas = deltas,
    merged_profiles = merged_profiles,
    final_input = working_input,
    processed_nodes = length(internal_nodes),
    chr = chr,
    allele = allele,
    branch_costs = branch_costs,
    new_edge_lengths = new_edge_lengths,
    history = history
  )
}

# Helper function to get post-order traversal of internal nodes
# get_post_order_internal_nodes <- function(tree) {
#   Ntip <- length(tree$tip.label)
#   Nnode <- tree$Nnode
#   all_internal_nodes <- (Ntip + 1):(Ntip + Nnode)

#   # Build children list
#   children_list <- vector("list", max(tree$edge))
#   for (i in seq_len(nrow(tree$edge))) {
#     parent <- tree$edge[i, 1]
#     child <- tree$edge[i, 2]
#     children_list[[parent]] <- c(children_list[[parent]], child)
#   }

#   # Post-order traversal
#   visited <- rep(FALSE, max(tree$edge))
#   post_order <- c()

#   post_order_traverse <- function(node) {
#     if (visited[node]) return()

#     # Visit children first
#     children <- children_list[[node]]
#     if (!is.null(children)) {
#       for (child in children) {
#         if (child > Ntip) {
#           # Only traverse internal nodes
#           post_order_traverse(child)
#         }
#       }
#     }

#     # Then visit this node
#     visited[node] <<- TRUE
#     if (node > Ntip) {
#       post_order <<- c(post_order, node)
#     }
#   }

#   # Start from root
#   root <- Ntip + 1
#   post_order_traverse(root)

#   post_order
# }

# Helper function to get direct children of a node
get_node_children <- function(tree, node) {
  children <- c()
  for (i in seq_len(nrow(tree$edge))) {
    if (tree$edge[i, 1] == node) {
      children <- c(children, tree$edge[i, 2])
    }
  }
  children
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

plot_tree_colored_by_bfb <- function(
  fit,
  chr,
  allele,
  B_dist,
  G_dist,
  add_legend = TRUE,
  add_title = TRUE,
  node_size = 3,
  tip_size = 2,
  branch_length = 1,
  colors = c("mediumpurple", "goldenrod"),
  show_tip_labels = FALSE
) {
  # Prepare tree
  tree <- fit$tree

  if (!ape::is.rooted(tree)) {
    tree <- phangorn::midpoint(tree)
  } else {
    tree <- phangorn::midpoint(tree)
  }

  # Get reconstruction results
  reconstruction <- reconstruct_tree(fit, chr, allele)

  # Extract and process delta values
  internal_nodes <- as.integer(names(reconstruction$deltas))
  delta_values <- unlist(reconstruction$deltas)

  # Create comprehensive node data
  node_data <- data.frame(
    node = internal_nodes,
    delta = delta_values,
    norm_delta = ifelse(delta_values > 0, 1, -1),
    decision = ifelse(delta_values > 0, "Likely BFB", "Other"),
    abs_delta = abs(delta_values),
    stringsAsFactors = FALSE,
    isTip = FALSE
  )

  # Handle zero-length branches
  if (!is.null(tree$edge.length)) {
    min_nonzero <- min(tree$edge.length[tree$edge.length > 0], na.rm = TRUE)
    if (!is.null(branch_length)) {
      tree$edge.length[tree$edge.length != branch_length] <- branch_length
    } else {
      tree$edge.length[tree$edge.length == 0] <- min_nonzero * 0.01
    }
  }

  # Create base tree plot
  p <- ggtree::ggtree(tree) %<+% node_data
  p <- p +
    ggplot2::theme(
      plot.margin = ggplot2::margin(20, 20, 20, 20),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )

  # Add tip points and labels
  if (show_tip_labels) {
    p <- p + ggtree::geom_tiplab(size = 3, offset = 0.01)
  }
  p <- p + ggtree::geom_tippoint(size = tip_size, color = "black", alpha = 0.7)

  # Add colored internal nodes based on coloring scheme
  if (nrow(node_data) > 0) {
    # Discrete coloring by decision
    color_map <- stats::setNames(colors, c("Other", "Likely BFB"))

    p <- p +
      ggtree::geom_nodepoint(
        ggplot2::aes(color = .data$decision),
        size = node_size,
        alpha = 1
      ) +
      ggplot2::scale_color_manual(
        values = color_map,
        name = "CN event"
      )
  }

  # Configure legend
  if (add_legend && nrow(node_data) > 0) {
    p <- p +
      ggplot2::guides(
        color = ggplot2::guide_legend(
          title = "Reconstruction Decision",
          title.position = "top",
          title.hjust = 0.5,
          override.aes = list(size = 4)
        )
      ) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = ggplot2::element_text(size = 12, face = "bold"),
        legend.text = ggplot2::element_text(size = 10),
        legend.key = ggplot2::element_rect(fill = "white", color = NA)
      )
  } else if (!add_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  # Add title
  if (add_title) {
    title_text <- paste0("Chr ", chr, ", Allele ", allele)
    p <- p +
      ggplot2::ggtitle(title_text) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 14,
          face = "bold",
          hjust = 0.5,
          margin = ggplot2::margin(b = 20)
        )
      )
  }

  return(p)
}


# Compare True BFB and Assigned BFB

compare_assignment_wrt_true = function(
  true_tree,
  cell_history,
  fit,
  allele,
  chr,
  B_dist,
  G_dist,
  colors = c("mediumpurple", "goldenrod"),
  node_size = 1
) {
  color_map <- stats::setNames(colors, c("Other", "Likely BFB"))

  cell_history$node_type = lapply(cell_history$cell_id, function(c_id) {
    d <- cell_history %>% dplyr::filter(.data$parent_id == c_id)
    if (nrow(d) != 0) {
      if (any(grepl("bfb", d$cn_event))) {
        return("Likely BFB")
      }
    }
    "Other"
  }) %>%
    unlist()

  true_tree_plot <- ggtree::ggtree(true_tree) %<+% cell_history
  true_tree_plot <- true_tree_plot +
    ggtree::geom_nodepoint(
      ggplot2::aes(color = .data$node_type),
      size = node_size,
      alpha = 1
    ) +
    ggplot2::scale_color_manual(values = color_map, name = "CN event") +
    ggplot2::ggtitle("True BFB") +
    ggplot2::theme(legend.position = "bottom")

  fit$tree <- true_tree
  reconstruction <- reconstruct_tree(fit, chr, allele)
  internal_nodes <- as.integer(names(reconstruction$deltas))
  delta_values <- unlist(reconstruction$deltas)
  node_data <- data.frame(
    node = internal_nodes,
    delta = delta_values,
    norm_delta = ifelse(delta_values > 0, 1, -1),
    decision = ifelse(delta_values > 0, "Likely BFB", "Other"),
    abs_delta = abs(delta_values),
    stringsAsFactors = FALSE,
    isTip = FALSE
  )
  p_tree_assigned <- ggtree::ggtree(true_tree) %<+% node_data
  p_tree_assigned <- p_tree_assigned +
    ggtree::geom_nodepoint(
      ggplot2::aes(color = .data$decision),
      size = node_size,
      alpha = 1
    ) +
    ggplot2::scale_color_manual(values = color_map, name = "CN event") +
    ggplot2::ggtitle("Assigned BFB") +
    ggplot2::theme(legend.position = "bottom")

  p_compare <- true_tree_plot | p_tree_assigned

  df <- lapply(internal_nodes, function(node_id) {
    pred <- node_data %>%
      dplyr::filter(.data$node == node_id) %>%
      dplyr::pull(.data$decision)
    true_cell_name <- c(true_tree$tip.label, true_tree$node.label)[node_id]
    true <- cell_history %>%
      dplyr::filter(.data$cell_id == true_cell_name) %>%
      dplyr::pull(.data$node_type)
    dplyr::bind_rows(node_id = node_id, pred = pred, true = true)
  }) %>%
    do.call("bind_rows", .)

  list(p_compare = p_compare, df = df)
}


compute_reconstructions = function(fit, chromosomes = NULL, alleles = NULL) {
  if (is.null(chromosomes)) {
    chromosomes = names(fit$all_input_Xs)
  }
  if (is.null(alleles)) {
    alleles = names(fit$all_input_Xs[[chromosomes[1]]])
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

  list(
    reconstructions = reconstructions,
    edge.length = new_edges,
    history = whole_history
  )
}
