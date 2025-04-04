
positive_selection_function = function(positive_selection_rate, hotspot_copies) {
  if (is.nan(hotspot_copies)) return(0)
  positive_selection_rate * (hotspot_copies > 1)
}

negative_selection_function = function(negative_selection_rate, hotspot_copies) {
  if (is.nan(hotspot_copies)) return(0)
  negative_selection_rate * !(hotspot_copies > 1)
}
