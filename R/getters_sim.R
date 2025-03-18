
get_seq_length_stats = function(sim) {
  seq_lenghts = lapply(sim$final_cells, function(s) {
    get_seq_length(s)
  }) %>% unlist()

  list(
    mean = mean(seq_lenghts),
    median = stats::median(seq_lenghts),
    sd = stats::sd(seq_lenghts)
  )
}
