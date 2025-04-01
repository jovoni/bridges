
cells2mat = function(cells, L, order) {
  m = matrix(0, nrow = length(cells), ncol = L)
  idx = 1
  for (idx in 1:length(cells)) {
    tab_vec = table(seq2vec(cells[[idx]]))
    m[idx, as.numeric(names(tab_vec))] = as.numeric(tab_vec)
  }
  rownames(m) = names(cells)

  if (order) {
    one_copy_percentage = rowSums(m == 1) / L

    # Get the order of the values in the ordering vector
    row_order <- order(one_copy_percentage, decreasing = TRUE)

    # Reorder the matrix rows
    m <- m[row_order, , drop = FALSE]
  }
  m
}
