
is_hotspot_gained = function(cell, hotspot) {
  table_vec = table(bridges:::seq2vec(cell))
  hotspot = c(hotspot)
  flag = all(lapply(hotspot, function(h) {h %in% names(table_vec)}) %>% unlist())
  if (flag) {
    all(table_vec[hotspot] > 1)
  } else {
    FALSE
  }
}
