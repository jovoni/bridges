
detect_bfb = function(fit, B_dist, G_dist, threshold = .005) {
  alleles = names(fit$all_input_Xs[[1]])
  chromosomes = names(fit$all_input_Xs)

  pseudo_cells_test_df = lapply(alleles, function(allele) {
    lapply(chromosomes, function(chr) {
      pseudo_cell_bin_test(fit = fit, chr = chr, allele = allele, B_dist = B_dist, G_dist = G_dist, threshold = threshold)
    }) %>% do.call("bind_rows", .)
  }) %>% do.call("bind_rows", .)

  pseudo_cells_test_df
}

pseudo_cell_t_test = function(fit, chr, allele, B_dist, G_dist, mu = 0, sign = F) {
  r = reconstruct_tree(fit, chr, allele, B_dist, G_dist)

  delta_values <- unlist(r$deltas)
  if (sign) delta_values = sign(delta_values)

  if (length(unique(delta_values)) == 1) {
    estimate = mean(delta_values)
    dplyr::tibble(mean = estimate, chr = chr, allele = allele, p.value=NA, greater = estimate > mu)
  } else {
    test = stats::t.test(delta_values, alternative = "greater", mu = mu)
    dplyr::tibble(mean = test$estimate, chr = chr, allele = allele, p.value=test$p.value, greater = test$estimate > mu)
  }
}

pseudo_cell_bin_test = function(fit, chr, allele, B_dist, G_dist, threshold) {
  r = reconstruct_tree(fit, chr, allele, B_dist, G_dist)

  delta_values <- unlist(r$deltas)
  N = length(delta_values)
  k = sum(delta_values > 0)

  test = stats::binom.test(k, N, p = threshold, alternative = "greater")
  dplyr::tibble(mean = test$estimate, N_trials = N, N_successes = k, chr = chr, allele = allele, p.value=test$p.value, threshold=threshold)
}


reconstruct_tree <- function(fit, chr, allele, B_dist, G_dist) {
  # Set target value based on allele type
  target_val <- ifelse(allele == "CN", 2, 1)

  # Get input data for the specified chromosome and allele
  input <- fit$all_input_Xs[[chr]][[allele]]

  # Get tree structure
  tree <- fit$tree

  # if (!ape::is.rooted(tree)) {
  #   tree = phangorn::midpoint(tree)
  # }
  tree$edge.length = rep(1, length(tree$edge.length))

  # Get post-order traversal of internal nodes to ensure children are processed before parents
  post_order_nodes <- get_post_order_internal_nodes(tree)
  if (is.null(tree$node.label)) tree$node.label <- paste0("N", seq_len(tree$Nnode))

  # Initialize results storage
  internal_nodes <- list()
  deltas <- list()
  merged_profiles = list()

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
        internal_node_cell <- b_d$ancestor
        merged_profiles[[as.character(node_id)]] = list(left = b_d$left, right = b_d$right)
      } else {
        internal_node_cell <- g_d$ancestor
        merged_profiles[[as.character(node_id)]]
        merged_profiles[[as.character(node_id)]] = NULL
      }

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
      working_input <- working_input[!rownames(working_input) %in% c(left_cell_name, right_cell_name),,drop=F]
      #print(rownames(working_input))
    }
  }

  # Return results as a list
  return(list(
    internal_nodes = internal_nodes,
    deltas = deltas,
    merged_profiles = merged_profiles,
    final_input = working_input,
    processed_nodes = length(internal_nodes),
    chr = chr, allele = allele
  ))
}

# Helper function to get post-order traversal of internal nodes
get_post_order_internal_nodes <- function(tree) {
  Ntip <- length(tree$tip.label)
  Nnode <- tree$Nnode
  all_internal_nodes <- (Ntip + 1):(Ntip + Nnode)

  # Build children list
  children_list <- vector("list", max(tree$edge))
  for (i in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    children_list[[parent]] <- c(children_list[[parent]], child)
  }

  # Post-order traversal
  visited <- rep(FALSE, max(tree$edge))
  post_order <- c()

  post_order_traverse <- function(node) {
    if (visited[node]) return()

    # Visit children first
    children <- children_list[[node]]
    if (!is.null(children)) {
      for (child in children) {
        if (child > Ntip) {  # Only traverse internal nodes
          post_order_traverse(child)
        }
      }
    }

    # Then visit this node
    visited[node] <<- TRUE
    if (node > Ntip) {
      post_order <<- c(post_order, node)
    }
  }

  # Start from root
  root <- Ntip + 1
  post_order_traverse(root)

  return(post_order)
}

# Helper function to get direct children of a node
get_node_children <- function(tree, node) {
  children <- c()
  for (i in seq_len(nrow(tree$edge))) {
    if (tree$edge[i, 1] == node) {
      children <- c(children, tree$edge[i, 2])
    }
  }
  return(children)
}

# Helper function to get clade pairs (as provided)
get_clade_pairs <- function(tree) {
  # Total number of tips and nodes
  Ntip <- length(tree$tip.label)
  Nnode <- tree$Nnode
  all_nodes <- (Ntip + 1):(Ntip + Nnode)

  # Build a list to store children of each node
  children_list <- vector("list", max(tree$edge))
  for (i in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    children_list[[parent]] <- c(children_list[[parent]], child)
  }

  # Recursive function to get all descendant tips of a node
  get_descendant_tips <- function(node) {
    if (node <= Ntip) {
      return(tree$tip.label[node])
    } else {
      children <- children_list[[node]]
      unlist(lapply(children, get_descendant_tips))
    }
  }

  # Store pairs of clades
  clade_pairs <- list()

  for (node in all_nodes) {
    children <- children_list[[node]]
    if (length(children) == 2) {
      left_tips <- get_descendant_tips(children[1])
      right_tips <- get_descendant_tips(children[2])
      clade_pairs[[as.character(node)]] <- list(left = left_tips, right = right_tips)
    }
  }

  return(clade_pairs)
}

plot_tree_colored_by_bfb <- function(fit, chr, allele, B_dist, G_dist,
                                     add_legend = TRUE, add_title = TRUE,
                                     node_size = 3, tip_size = 2,
                                     branch_length = 1,
                                     colors = c("mediumpurple", "goldenrod"),
                                     show_tip_labels = FALSE) {

  # Prepare tree
  tree <- fit$tree

  if (!ape::is.rooted(tree)) {
    tree <- phangorn::midpoint(tree)
  } else {
    tree <- phangorn::midpoint(tree)
  }

  # Get reconstruction results
  reconstruction <- reconstruct_tree(fit, chr, allele, B_dist, G_dist)

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

compare_assignment_wrt_true = function(true_tree, cell_history,
                                       fit, allele, chr,
                                       B_dist, G_dist,
                                       colors = c("mediumpurple", "goldenrod"),
                                       node_size = 1) {

  color_map <- stats::setNames(colors, c("Other", "Likely BFB"))

  # True Tree plot
  cell_history$node_type = lapply(cell_history$cell_id, function(c_id) {
    d <- cell_history %>% dplyr::filter(.data$parent_id == c_id)
    if (nrow(d) != 0) {
      if ("bfb" %in% unique(d$cn_event) ) {
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
    ggplot2::scale_color_manual(
      values = color_map,
      name = "CN event"
    ) +
    ggplot2::ggtitle("True BFB") +
    ggplot2::theme(legend.position = "bottom")


  # Assigned Tree plot

  # Get reconstruction results
  fit$tree <- true_tree
  reconstruction <- reconstruct_tree(fit, chr, allele, B_dist, G_dist)

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

  # Create assigned tree
  p_tree_assigned <- ggtree::ggtree(true_tree) %<+% node_data
  p_tree_assigned <- p_tree_assigned +
    ggtree::geom_nodepoint(
      ggplot2::aes(color = .data$decision),
      size = node_size,
      alpha = 1
    ) +
    ggplot2::scale_color_manual(
      values = color_map,
      name = "CN event"
    ) +
    ggplot2::ggtitle("Assigned BFB") +
    ggplot2::theme(legend.position = "bottom")

  p_compare <- true_tree_plot | p_tree_assigned

  # Extract values
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
