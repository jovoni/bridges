
# Takes a CNA data format and returns set of breakpoints
get_breakpoints_dist = function(cna_matrix, allele) {
  lapply(1:nrow(cna_matrix), function(j) {
    which(diff(cna_matrix[j,]) != 0)
  }) %>% unlist()
}

get_cn_length_rate = function(cna_matrix, allele) {
  cna_lengths = lapply(1:nrow(cna_matrix), function(j) {
    idxs = which(diff(cna_matrix[j,]) != 0)
    if (length(idxs) > 1) {
      return(diff(idxs))
    } else {
      return(NULL)
    }
  }) %>% unlist()
  return(mean(cna_lengths))
}

get_gain_profile = function(cna_matrix, base_value) {
  colMeans(cna_matrix > base_value)
}

get_loss_profile = function(cna_matrix, base_value) {
  colMeans(cna_matrix < base_value)
}

get_avg_cn_profile = function(cna_matrix) {
  colMeans(cna_matrix)
}

run_ABC = function(cna_data, allele, chromosome, chr, pos, params) {
  # Extract parameters
  amp_rate = params$amp_rate
  del_rate = params$del_rate
  bfb_prob = params$bfb_prob
  birth_rate = params$birth_rate
  death_rate = params$death_rate
  lambda = params$lambda

  base_value = ifelse(allele == "CN", 2, 1)
  target_cna_matrix = bridges:::tibble_to_matrix(cna_data, chromosome = chromosome, value_column = allele)
  target_gain_dist = get_gain_profile(target_cna_matrix, base_value = base_value)
  target_loss_dist = get_loss_profile(target_cna_matrix, base_value = base_value)
  target_avg_cn_dist = get_avg_cn_profile(target_cna_matrix)
  n_target_cells = nrow(target_cna_matrix)
  target_cna_rate = get_cn_length_rate(target_cna_matrix, allele)
  custom_breakpoint_dist = get_breakpoints_dist(target_cna_matrix, allele = allele)

  # Run simulation with provided parameters
  sim = bridge_sim(initial_cells = 1,
                   bin_length = 5e5,
                   allow_wgd = TRUE,
                   max_time = 300,
                   max_cells = n_target_cells,
                   normal_dup_rate = 0,
                   amp_rate = amp_rate,        # Now using parameter
                   del_rate = del_rate,        # Now using parameter
                   bfb_prob = bfb_prob,        # Now using parameter
                   birth_rate = birth_rate,    # Now using parameter
                   death_rate = death_rate,    # Now using parameter
                   rate = target_cna_rate,
                   lambda = lambda,            # Now using parameter
                   bfb_allele = paste0(chromosome, ":", allele),
                   chromosomes = c(chr),
                   first_round_of_bfb = TRUE,
                   positive_selection_rate = 0,
                   negative_selection_rate = 0,
                   hotspot = list(chr = paste0(chr, ":", allele), pos = pos),
                   breakpoint_support = "custom",
                   custom_breakpoints = custom_breakpoint_dist)

  # Calculate summary statistics
  sim_cna_mat = bridges:::tibble_to_matrix(sim$cna_data, chromosome = chromosome, value_column = allele)
  pred_gain_dist = get_gain_profile(sim_cna_mat, base_value = base_value)
  pred_loss_dist = get_loss_profile(sim_cna_mat, base_value = base_value)
  pred_avg_cn_dist = get_avg_cn_profile(sim_cna_mat)

  # Calculate and return distance metrics
  rmse_gain = caret::RMSE(pred_gain_dist, target_gain_dist)
  rmse_loss = caret::RMSE(pred_loss_dist, target_loss_dist)
  rmse_avg_cn = caret::RMSE(pred_avg_cn_dist, target_avg_cn_dist)

  # Combined distance metric (weighted sum)
  total_distance = rmse_gain + rmse_loss + rmse_avg_cn

  return(list(
    distance = total_distance,
    rmse_gain = rmse_gain,
    rmse_loss = rmse_loss,
    rmse_avg_cn = rmse_avg_cn,
    params = params
  ))
}

# Prior sampling function
sample_priors = function() {
  # Sample relative proportions for amp_rate, del_rate, bfb_prob from Dirichlet
  # Using symmetric Dirichlet with alpha = 1 (uniform on simplex)
  # You can adjust alpha values to encode prior beliefs about relative rates
  dirichlet_alpha = c(1, 1, 1)  # Equal prior weight for amp, del, bfb
  relative_rates = MCMCpack::rdirichlet(1, dirichlet_alpha)[1, ]

  # Use relative rates directly (they sum to 1)
  amp_rate = relative_rates[1]
  del_rate = relative_rates[2]
  bfb_prob = relative_rates[3]

  # Fix birth rate to 1 and sample death rate as fraction of birth rate
  birth_rate = 1.0
  death_fraction = runif(1, 0, 0.8)  # Death rate as fraction of birth rate [0, 0.8]
  death_rate = birth_rate * death_fraction

  list(
    amp_rate = amp_rate,
    del_rate = del_rate,
    bfb_prob = bfb_prob,
    birth_rate = birth_rate,
    death_rate = death_rate,
    lambda = runif(1, 0.001, 0.1),      # Keep lambda as before
    death_fraction = death_fraction     # Store for analysis
  )
}

# Main ABC algorithm
abc_inference = function(cna_data, allele, chromosome, chr, pos,
                         n_simulations = 10000,
                         tolerance_quantile = 0.01,
                         n_cores = 1) {

  # Storage for results
  results = vector("list", n_simulations)

  # Progress bar
  pb = utils::txtProgressBar(min = 0, max = n_simulations, style = 3)

  # Run simulations
  if (n_cores > 1) {

    results = parallel::mclapply(1:n_simulations, function(i) {
      params = sample_priors()
      tryCatch({
        run_ABC(cna_data, allele, chromosome, chr, pos, params)
      }, error = function(e) {
        list(distance = Inf, params = params, error = e$message)
      })
    }, mc.cores = n_cores)

  } else {
    # Sequential execution
    for (i in 1:n_simulations) {
      params = sample_priors()

      results[[i]] = tryCatch({
        run_ABC(cna_data, allele, chromosome, chr, pos, params)
      }, error = function(e) {
        list(distance = Inf, params = params, error = e$message)
      })


      utils::setTxtProgressBar(pb, i)
    }
  }
  close(pb)

  # Filter out failed simulations
  valid_results = results[sapply(results, function(x) is.finite(x$distance))]

  if (length(valid_results) == 0) {
    stop("No valid simulations completed. Check your model parameters and data.")
  }

  # Extract distances and sort
  distances = sapply(valid_results, function(x) x$distance)
  sorted_indices = order(distances)

  # Apply tolerance threshold
  n_accept = max(1, floor(length(valid_results) * tolerance_quantile))
  accepted_indices = sorted_indices[1:n_accept]
  accepted_results = valid_results[accepted_indices]

  # Extract accepted parameters
  accepted_params = data.frame(
    amp_rate = sapply(accepted_results, function(x) x$params$amp_rate),
    del_rate = sapply(accepted_results, function(x) x$params$del_rate),
    bfb_prob = sapply(accepted_results, function(x) x$params$bfb_prob),
    birth_rate = sapply(accepted_results, function(x) x$params$birth_rate),
    death_rate = sapply(accepted_results, function(x) x$params$death_rate),
    lambda = sapply(accepted_results, function(x) x$params$lambda),
    death_fraction = sapply(accepted_results, function(x) x$params$death_fraction),
    distance = sapply(accepted_results, function(x) x$distance)
  )

  # Calculate summary statistics
  param_summary = data.frame(
    parameter = names(accepted_params)[-8], # Exclude distance column
    mean = sapply(accepted_params[,-8], stats::mean),
    median = sapply(accepted_params[,-8], stats::median),
    sd = sapply(accepted_params[,-8], stats::sd),
    q025 = sapply(accepted_params[,-8], stats::quantile, 0.025),
    q975 = sapply(accepted_params[,-8], stats::quantile, 0.975)
  )

  return(list(
    accepted_params = accepted_params,
    param_summary = param_summary,
    n_simulations = n_simulations,
    n_accepted = n_accept,
    acceptance_rate = n_accept / length(valid_results),
    tolerance_threshold = max(accepted_params$distance)
  ))
}

# Diagnostic plotting function
plot_abc_results = function(abc_results) {
  library(ggplot2)
  library(gridExtra)

  params = abc_results$accepted_params

  plots = list()

  # Plot the three relative rates (they sum to 1)
  relative_params = c("amp_rate", "del_rate", "bfb_prob")
  for (i in 1:3) {
    param_name = relative_params[i]


    p = ggplot2::ggplot(params, ggplot2::aes_string(x = param_name)) +
      ggplot2::geom_histogram(bins = 30, alpha = 0.7, fill = "steelblue") +
      ggplot2::geom_vline(xintercept = stats::median(params[[param_name]]),
                 color = "red", linetype = "dashed") +
      ggplot2::labs(title = paste("Relative:", param_name),
           x = "Relative frequency", y = "Count") +
      ggplot2::xlim(0, 1) +
      ggplot2::theme_minimal()

    plots[[i]] = p
  }

  # Plot other parameters
  other_params = c("death_rate", "death_fraction", "lambda")
  for (i in 1:3) {
    param_name = other_params[i]

    p = ggplot2::ggplot(params, ggplot2::aes_string(x = param_name)) +
      ggplot2::geom_histogram(bins = 30, alpha = 0.7, fill = "darkgreen") +
      ggplot2::geom_vline(xintercept = stats::median(params[[param_name]]),
                 color = "red", linetype = "dashed") +
      ggplot2::labs(title = paste("Posterior:", param_name),
           x = param_name, y = "Count") +
      ggplot2::theme_minimal()

    plots[[i + 3]] = p
  }

  # Ternary-like plot showing the three relative rates
  # Since they sum to 1, we can show relationships between pairs
  p = ggplot2::ggplot(params, ggplot2::aes(x = amp_rate, y = del_rate)) +
    ggplot2::geom_point(alpha = 0.6, ggplot2::aes(color = bfb_prob)) +
    ggplot2::scale_color_gradient(low = "blue", high = "red", name = "BFB prob") +
    ggplot2::labs(title = "Relative Rates (Amp vs Del, colored by BFB)",
         x = "Amp rate", y = "Del rate") +
    ggplot2::theme_minimal()
  plots[[7]] = p

  # Death fraction vs relative rates
  p = ggplot2::ggplot(params, ggplot2::aes(x = death_fraction, y = amp_rate)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::labs(title = "Death Fraction vs Amp Rate",
         x = "Death Fraction", y = "Amp Rate") +
    ggplot2::theme_minimal()
  plots[[8]] = p

  # Distance vs main parameters
  p = ggplot2::ggplot(params, ggplot2::aes(x = amp_rate, y = distance)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "red") +
    ggplot2::labs(title = "Distance vs Amp Rate", x = "Amp Rate", y = "Distance") +
    ggplot2::theme_minimal()
  plots[[9]] = p

  gridExtra::grid.arrange(grobs = plots, ncol = 3)
}
