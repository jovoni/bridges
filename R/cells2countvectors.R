
cells2countvectors = function(cells, L) {
  breakpoints = get_breakpoints(cells, L)
  breakpoints = c(0, breakpoints, L)

  list_countvectors = lapply(cells, function(cell) {
    count_vector = lapply(2:length(breakpoints), function(i) {
      b_prev = breakpoints[i-1]
      b_curr = breakpoints[i]

      vec = seq2vec(cell)
      tab_vec = table(vec[vec > b_prev & vec <= b_curr])
      if (length(tab_vec) == 0) {
        count = 0
      } else {
        count = unique(tab_vec)
      }

      if (length(count) > 1) {
        stop("Error with current breakpoints.")
      }

      count
    }) %>% unlist()
  })

  do.call("rbind", list_countvectors)
}
