
compute_distance_matrix_enhanced <- function(A, B = NULL, alpha = NULL, cn_weight = NULL) {
  if (is.null(B)) {
    if (is.vector(A)) {
      A <- matrix(A, nrow = 1)
    }

    # Compute pairwise distances within A
    N <- nrow(A)
    D <- matrix(0, nrow = N, ncol = N)
    if (N == 1) return(D)

    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        D[i, j] <- custom_distance(A[i, ], A[j, ], alpha = alpha, cn_weight = cn_weight)
        D[j, i] <- D[i, j]  # Symmetric matrix
      }
    }
    return(D)
  } else {
    if (is.vector(A)) {
      A <- matrix(A, nrow = 1)
    }
    if (is.vector(B)) {
      B <- matrix(B, nrow = 1)
    }

    # Compute distances between each row of A and each row of B
    N <- nrow(A)
    M <- nrow(B)
    D <- matrix(0, nrow = N, ncol = M)
    for (i in 1:N) {
      for (j in 1:M) {
        D[i, j] <- custom_distance(A[i, ], B[j, ], alpha = alpha, cn_weight = cn_weight)
      }
    }
    return(D)
  }
}

propose_groups_of_cells = function(D, min_dist) {
  pairs <- which(D == min_dist & lower.tri(D), arr.ind = TRUE)
  pairs <- cbind(rownames(D)[pairs[,1]], colnames(D)[pairs[,2]])

  # Convert to list of named sets
  pair_sets <- apply(pairs, 1, function(x) list(sort(x)))
  pair_sets <- lapply(pair_sets, function(x) x[[1]])

  return(pair_sets)

  # Merge overlapping sets
  if (length(pair_sets) > 0) {
    changed <- TRUE
    while (changed) {
      changed <- FALSE
      new_sets <- list()

      for (current in pair_sets) {
        merged <- FALSE
        for (i in seq_along(new_sets)) {
          if (length(intersect(current, new_sets[[i]]))) {
            new_sets[[i]] <- union(new_sets[[i]], current)
            merged <- changed <- TRUE
            break
          }
        }
        if (!merged) new_sets <- c(new_sets, list(current))
      }
      pair_sets <- new_sets
    }
  }

  # Obtain only maximal vector in each pair set
  otus = lapply(pair_sets, extract_maximal_cells)

  lapply(otus, sort)
}

clean_group_of_cells = function(otus, C) {
  clean_otu = function(otu) {
    lapply(otu, function(o) {
      s = sum(grepl(gsub("\\|", "-", o), gsub("\\|", "-", otu), fixed = TRUE))
      if (s == 1) {return(o)}
    }) %>% unlist()
  }

  otus = lapply(otus, clean_otu)

  otus = lapply(otus, function(otu) {
    if (length(otu) == 1) return(NULL)
    #if (all(otu %in% unlist(cells_seen_together))) return(NULL)
    #otu_name = paste0(otu, collapse = "")
    #otu_name = paste0("|", otu_name, "|")

    #f = otu_name %in% unlist(cells_seen_together)

    return(otu)
  })

  otus[sapply(otus, is.null)] <- NULL

  # if (length(otus)) {
  #   new_pseudocells = lapply(otus, function(nu) {
  #     build_new_pseudocell(nu, C)
  #   }) %>% do.call("rbind", .)
  #
  #
  #   kept_indices <- which(!duplicated(new_pseudocells))
  #   new_pseudocells = new_pseudocells %>% unique()
  #   otus = lapply(kept_indices, function(i) otus[[i]])
  #
  #   # Check that there is no repeated cell in good_ones and if produce a integer representaiton
  #   good_idxs = which(lapply(1:nrow(new_pseudocells), function(idx) {
  #     ps = rownames(new_pseudocells)[idx]
  #     cells_used = unlist(strsplit(gsub("\\|", "", ps), "-"))
  #     y = new_pseudocells[idx,]
  #     (length(cells_used) == length(unique(cells_used)) & !(ps %in% rownames(C))) & all(y == floor(y))
  #   }) %>% unlist() == TRUE)
  #
  #   #if (length(good_idxs_2) == 0) {print(pippo)}
  #
  #   good_ones = rownames(new_pseudocells)[good_idxs]
  #
  #   if (length(good_ones)) {
  #     otus = lapply(good_idxs, function(idx) {otus[[idx]]})
  #     print("FOUND")
  #     parents = matrix(new_pseudocells[good_idxs,], nrow = length(good_ones), ncol = ncol(new_pseudocells))
  #     rownames(parents) = good_ones
  #     return(list(daughters = otus, parents = parents))
  #   }
  # }

  if (length(otus)) {
    # Build new pseudocells
    new_pseudocells = lapply(otus, function(nu) {
      build_new_pseudocell(nu, C)
    }) %>% do.call("rbind", .)

    # Ensure new_pseudocells is always a matrix, even if it has only one row
    if (!is.matrix(new_pseudocells)) {
      new_pseudocells <- matrix(new_pseudocells, nrow = 1)
      rownames(new_pseudocells) <- names(otus)
    }

    # First filter: check for bad pseudo_cells
    good_idxs = which(lapply(1:nrow(new_pseudocells), function(idx) {
      ps = rownames(new_pseudocells)[idx]
      cells_used = unlist(strsplit(gsub("\\|", "", ps), "-"))
      y = new_pseudocells[idx,]
      (length(cells_used) == length(unique(cells_used)) & !(ps %in% rownames(C))) & all(y == floor(y))
    }) %>% unlist() == TRUE)

    # Only proceed if we have good pseudo_cells
    if (length(good_idxs) > 0) {
      # Keep only the good pseudocells and corresponding otus
      new_pseudocells = new_pseudocells[good_idxs, , drop = FALSE]  # drop=FALSE preserves matrix structure
      otus = lapply(good_idxs, function(idx) otus[[idx]])

      # Now remove duplicates from the filtered pseudocells
      kept_indices <- which(!duplicated(new_pseudocells))
      new_pseudocells = new_pseudocells[kept_indices, , drop = FALSE]  # drop=FALSE preserves matrix structure

      # Keep only the otus that correspond to the distinct pseudocells
      otus = lapply(kept_indices, function(i) otus[[i]])

      good_ones = rownames(new_pseudocells)
      if (length(good_ones)) {
        #print("FOUND")
        parents = matrix(new_pseudocells, nrow = length(good_ones), ncol = ncol(new_pseudocells))
        rownames(parents) = good_ones
        return(list(daughters = otus, parents = parents))
      }
    }
  }
}


extend_M = function(M, C) {
  new_N = nrow(C)
  old_N = nrow(M)
  M_extended <- matrix(FALSE, nrow = new_N, ncol = new_N)  # Initialize new matrix with FALSE
  M_extended[1:old_N, 1:old_N] <- M  # Copy original values
  M = M_extended
  diag(M) = TRUE
  M
}


get_progeny <- function(tree, parent) {
  if (!parent %in% names(tree)) {
    return(character(0))  # Return empty if parent not found in the tree
  }

  direct_children <- tree[[parent]]  # Get direct children
  all_descendants <- direct_children  # Initialize with direct children

  for (child in direct_children) {
    all_descendants <- c(all_descendants, get_progeny(tree, child))  # Recursively add progeny
  }

  return(unique(all_descendants))  # Return unique values to avoid cycles
}

update_M = function(M, daughters, pd_list, all = TRUE) {

  for (i in 1:length(daughters)) {
    parent = names(daughters)[i]
    daughter = daughters[[i]]

    progeny = get_progeny(pd_list, parent)
    progeny_and_parent = c(progeny, parent)

    # Update the matrix to set TRUE for all pairwise combinations
    idxs = which(rownames(M) %in% progeny_and_parent)
    M[idxs, idxs] <- TRUE
  }

  M
}

