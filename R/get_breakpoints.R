
# Function to obtain cuts in the profile
get_breakpoints = function(cells, L) {
  m = cells2mat(cells, L, order = FALSE)

  breakpoints = lapply(2:ncol(m), function(i) {
    if (any(m[,i] != m[,i-1])) {
      return(i-1)
    }
  }) %>% unlist()
  breakpoints
}
