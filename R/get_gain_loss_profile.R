
get_gain_loss_profile = function(cells, L) {
  m = matrix(0, nrow = length(cells), ncol = L)
  for (idx in 1:length(cells)) {
    tab_vec = table(seq2vec(cells[[idx]]))
    m[idx, as.numeric(names(tab_vec))] = as.numeric(tab_vec)
  }

  list(Loss = -colMeans(m < 1), Gain = colMeans(m * (m > 1)))
}
