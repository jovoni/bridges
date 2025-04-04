
is_hotspot_gained = function(cell, hotspot) {
  hc = get_hotspot_copies(cell, hotspot)
  if (is.nan(hc)) return(NA)
  hc > 1
}

get_hotspot_copies = function(cell, hotspot) {
  if (is.null(hotspot)) return(NaN)

  table_vec = table(seq2vec(cell))
  hotspot = c(hotspot)
  flag = names(table_vec) %in% hotspot
  if (any(flag)) {
    table_vec[names(table_vec) %in% hotspot] %>% mean()
  } else {
    0
  }
}
