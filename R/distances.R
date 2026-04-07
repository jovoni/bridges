# G distances ####
# Original one
find_greedy_distance = function(a, b, target_val = 2) {
  # Make copies to avoid modifying inputs
  current = a
  x = current - b

  # Initialize history with starting vector
  history = list(current)
  cost = 0

  while (!all(x == 0)) {
    # Find first position where vectors differ
    if (any(x != 0)) {
      head = min(which(x != 0))

      # Find the end of the contiguous block with same sign
      end = head
      sign_val = sign(x[head])

      while (end < length(x) && sign(x[end + 1]) == sign_val) {
        end = end + 1
      }

      #step = min(abs(x[head:end]))
      # Apply the operation (add or subtract 1)
      if (x[head] < 0) {
        # Need to increase current to match b
        #current[head:end] = current[head:end] + step
        current[head:end] = current[head:end] + 1
      } else {
        # Need to decrease current to match b
        #current[head:end] = current[head:end] - step
        current[head:end] = current[head:end] - 1
      }

      # Update history and cost
      history = c(history, list(current))
      cost = cost + 1

      # Recalculate difference
      x = current - b
    }
  }

  # Calculate ancestor - vector closest to target_val
  distances_to_target = sapply(history, function(x) {
    sum(abs(x - target_val))
  })
  ancestor_index = which.min(distances_to_target)
  ancestor = history[[ancestor_index]]
  b_cost = length(history) - ancestor_index
  a_cost = cost - b_cost

  # starting_vec_list = list(a=a, b=b)
  # distances_to_target = sapply(starting_vec_list, function(x) {sum(abs(x - target_val))})
  # ancestor = starting_vec_list[[which.min(distances_to_target)]]

  return(list(
    history = history,
    cost = cost,
    ancestor = ancestor,
    b_cost = b_cost,
    a_cost = a_cost,
    a_delta = a - ancestor,
    b_delta = b - ancestor
  ))
}

find_greedy_distance_with_steps = function(a, b, target_val = 2) {
  # Make copies to avoid modifying inputs
  current = a
  x = current - b

  # Initialize history with starting vector
  history = list(current)
  cost = 0

  while (!all(x == 0)) {
    # Find first position where vectors differ
    if (any(x != 0)) {
      head = min(which(x != 0))

      # Find the end of the contiguous block with same sign
      end = head
      sign_val = sign(x[head])

      while (end < length(x) && sign(x[end + 1]) == sign_val) {
        end = end + 1
      }

      step = min(abs(x[head:end]))
      # Apply the operation (add or subtract 1)
      if (x[head] < 0) {
        # Need to increase current to match b
        current[head:end] = current[head:end] + step
        #current[head:end] = current[head:end] + 1
      } else {
        # Need to decrease current to match b
        current[head:end] = current[head:end] - step
        #current[head:end] = current[head:end] - 1
      }

      # Update history and cost
      history = c(history, list(current))
      cost = cost + 1

      # Recalculate difference
      x = current - b
    }
  }

  # Calculate ancestor - vector closest to target_val

  distances_to_target = sapply(history, function(x) sum(abs(x - target_val)))
  ancestor_index = which.min(distances_to_target)
  ancestor = history[[ancestor_index]]
  b_cost = length(history) - ancestor_index
  a_cost = cost - b_cost

  # starting_vec_list = list(a=a, b=b)
  # distances_to_target = sapply(starting_vec_list, function(x) {sum(abs(x - target_val))})
  # ancestor = starting_vec_list[[which.min(distances_to_target)]]

  return(list(
    history = history,
    cost = cost,
    ancestor = ancestor,
    b_cost = b_cost,
    a_cost = a_cost,
    a_delta = a - ancestor,
    b_delta = b - ancestor
  ))
}

# B matrices ####

find_greedy_distance_with_avg = function(a, b, penalty = 0) {
  # Make a copy of input vector to avoid modifying the original
  current = a

  # Calculate initial difference and set even differences to 0
  x = current - b
  x[x %% 2 == 0] = 0

  # Initialize history and cost
  history = list(current)
  cost = 0

  # Continue until all differences are even (or zero)
  while (!all(x %% 2 == 0)) {
    # Find first position with non-zero difference
    if (any(x != 0)) {
      head = min(which(x != 0))

      # Find the end of the contiguous block with same sign
      end = head
      sign_val = sign(x[head])

      # Make sure not to go out of bounds
      while (end < length(x) && sign(x[end + 1]) == sign_val) {
        end = end + 1
      }

      # Apply the operation (add or subtract 1)
      if (x[head] < 0) {
        current[head:end] = current[head:end] + 1
      } else {
        current[head:end] = current[head:end] - 1
      }

      # Update history and cost
      history = c(history, list(current))
      cost = cost + 1

      # Recalculate difference and set even differences to 0
      x = current - b
      x[x %% 2 == 0] = 0
    }
  }

  # Calculate average of final state and target
  avg = (current + b) / 2
  a_cost = cost + 0.5 + penalty / 2
  b_cost = 0.5 + penalty / 2

  # Add average and final state to history
  history = c(history, list(avg), list(b))
  cost = cost + 1 + penalty

  return(list(
    history = history,
    cost = cost,
    ancestor = avg,
    b_cost = b_cost,
    a_cost = a_cost,
    a_delta = a - avg,
    b_delta = b - avg,
    left = current,
    right = b
  ))
}

find_greedy_distance_A = function(a, b, penalty = 0) {
  # Make a copy of input vector to avoid modifying the original
  current = a

  # Calculate initial difference
  x = current - b

  # Initialize history and cost
  history = list(current)
  cost = 0

  # Continue until all differences are even (or zero)
  while (!all(x %% 2 == 0)) {
    # Find positions with odd differences
    odd_positions = which(x %% 2 != 0)

    if (length(odd_positions) > 0) {
      # Find first position with odd difference
      head = min(odd_positions)

      # Find the end of the contiguous block with same sign among odd positions
      end = head
      sign_val = sign(x[head])

      # Extend to consecutive positions with same sign (considering all positions, not just odd)
      while (end < length(x) && sign(x[end + 1]) == sign_val) {
        end = end + 1
      }

      # Apply the operation (add or subtract 1) to the entire block
      if (x[head] < 0) {
        current[head:end] = current[head:end] + 1
      } else {
        current[head:end] = current[head:end] - 1
      }

      # Update history and cost
      history = c(history, list(current))
      cost = cost + 1

      # Recalculate difference (keep all values, don't zero out even ones)
      x = current - b
    }
  }

  # Calculate average of final state and target
  avg = (current + b) / 2
  a_cost = cost + 0.5 + penalty / 2
  b_cost = 0.5 + penalty / 2

  # Add average and final state to history
  history = c(history, list(avg), list(b))
  cost = cost + 1 + penalty

  return(list(
    history = history,
    cost = cost,
    ancestor = avg,
    b_cost = b_cost,
    a_cost = a_cost,
    a_delta = a - avg,
    b_delta = b - avg,
    left = current,
    right = b
  ))
}

find_greedy_distance_A_contig = function(a, b, penalty = 0) {
  # Make a copy of input vector to avoid modifying the original
  current = a

  # Calculate initial difference
  x = current - b

  # Initialize history and cost
  history = list(current)
  cost = 0

  # Continue until all differences are even (or zero)
  while (!all(x %% 2 == 0)) {
    #print(x)

    # Find positions with odd differences
    odd_positions = which(x %% 2 != 0)

    if (length(odd_positions) > 0) {
      # Find first position with odd difference
      head = min(odd_positions)

      # Find the end of the contiguous block with same sign among odd positions
      end = head
      sign_val = sign(x[head])

      # Extend to consecutive positions with same sign (considering all positions, not just odd)
      while (end < length(x) && sign(x[end + 1]) == sign_val) {
        end = end + 1
      }

      # Apply the operation (add or subtract 1) to the entire block
      if (x[head] < 0) {
        current[head:end] = current[head:end] + 1
      } else {
        current[head:end] = current[head:end] - 1
      }

      # Update history and cost
      history = c(history, list(current))
      cost = cost + 1

      # Recalculate difference (keep all values, don't zero out even ones)
      x = current - b
    }
  }

  # Calculate average of final state and target
  r = rle(x != 0)
  avg = (current + b) / 2
  cost = cost + sum(r$values) + penalty

  # if (all(x == 0)) cost = Inf

  avg = (current + b) / 2
  a_cost = cost - 1
  b_cost = 1

  # Add average and final state to history
  history = c(history, list(avg), list(b))

  return(list(
    history = history,
    cost = cost,
    ancestor = avg,
    b_cost = b_cost,
    a_cost = a_cost,
    a_delta = a - avg,
    b_delta = b - avg,
    left = current,
    right = b
  ))
}

# find_greedy_distance_A_each = function(a, b, penalty = 0) {
#   # Make a copy of input vector to avoid modifying the original
#   current = a
#
#   # Calculate initial difference
#   x = current - b
#
#   # Initialize history and cost
#   history = list(current)
#   cost = 0
#
#   # Continue until all differences are even (or zero)
#   while(!all(x %% 2 == 0)) {
#     #print(x)
#
#     # Find positions with odd differences
#     odd_positions = which(x %% 2 != 0)
#
#     if(length(odd_positions) > 0) {
#       # Find first position with odd difference
#       head = min(odd_positions)
#
#       # Find the end of the contiguous block with same sign among odd positions
#       end = head
#       sign_val = sign(x[head])
#
#       # Extend to consecutive positions with same sign (considering all positions, not just odd)
#       while(end < length(x) && sign(x[end + 1]) == sign_val) {
#         end = end + 1
#       }
#
#       # Apply the operation (add or subtract 1) to the entire block
#       if(x[head] < 0) {
#         current[head:end] = current[head:end] + 1
#       } else {
#         current[head:end] = current[head:end] - 1
#       }
#
#       # Update history and cost
#       history = c(history, list(current))
#       cost = cost + 1
#
#       # Recalculate difference (keep all values, don't zero out even ones)
#       x = current - b
#     }
#   }
#
#   # Calculate average of final state and target
#   avg = (current + b) / 2
#   cost = cost + sum(x != 0) + penalty
#
#   # Add average and final state to history
#   history = c(history, list(avg), list(b))
#
#   # Return results (assuming the original function returns something)
#   return(list(
#     final_state = b,
#     cost = cost,
#     history = history,
#     ancestor = avg
#   ))
# }

# find_greedy_distance_B = function(a, b, penalty = 0) {
#   # Make a copy of input vector to avoid modifying the original
#   current = a
#
#   # Calculate initial difference
#   x = current - b
#
#   # Initialize history and cost
#   history = list(current)
#   cost = 0
#
#   # Continue until all differences are even (or zero)
#   while(!all(x %% 2 == 0)) {
#     #print(x)
#
#     # Find positions with odd differences
#     odd_positions = which(x %% 2 != 0)
#
#     if(length(odd_positions) > 0) {
#       # Find first position with odd difference
#       head = min(odd_positions)
#
#       # Find the end of the contiguous block (ignoring sign)
#       end = head
#
#       # Extend to consecutive positions with non-zero differences (regardless of sign)
#       while(end < length(x) && x[end + 1] != 0) {
#         end = end + 1
#       }
#
#       # Apply the operation based on the sign at the head position
#       if(x[head] < 0) {
#         current[head:end] = current[head:end] + 1
#       } else {
#         current[head:end] = current[head:end] - 1
#       }
#
#       # Update history and cost
#       history = c(history, list(current))
#       cost = cost + 1
#
#       # Recalculate difference (keep all values, don't zero out even ones)
#       x = current - b
#     }
#   }
#
#   # Calculate average of final state and target
#   avg = (current + b) / 2
#   cost = cost + 1 + penalty
#
#   # Add average and final state to history
#   history = c(history, list(avg), list(b))
#
#   # Return results (assuming the original function returns something)
#   return(list(
#     final_state = b,
#     cost = cost,
#     history = history,
#     ancestor = avg
#   ))
# }

# find_greedy_distance_B_each = function(a, b, penalty = 0) {
#   # Make a copy of input vector to avoid modifying the original
#   current = a
#
#   # Calculate initial difference
#   x = current - b
#
#   # Initialize history and cost
#   history = list(current)
#   cost = 0
#
#   # Continue until all differences are even (or zero)
#   while(!all(x %% 2 == 0)) {
#     #print(x)
#
#     # Find positions with odd differences
#     odd_positions = which(x %% 2 != 0)
#
#     if(length(odd_positions) > 0) {
#       # Find first position with odd difference
#       head = min(odd_positions)
#
#       # Find the end of the contiguous block (ignoring sign)
#       end = head
#
#       # Extend to consecutive positions with non-zero differences (regardless of sign)
#       while(end < length(x) && x[end + 1] != 0) {
#         end = end + 1
#       }
#
#       # Apply the operation based on the sign at the head position
#       if(x[head] < 0) {
#         current[head:end] = current[head:end] + 1
#       } else {
#         current[head:end] = current[head:end] - 1
#       }
#
#       # Update history and cost
#       history = c(history, list(current))
#       cost = cost + 1
#
#       # Recalculate difference (keep all values, don't zero out even ones)
#       x = current - b
#     }
#   }
#
#   # Calculate average of final state and target
#   avg = (current + b) / 2
#   #cost = cost + 1 + penalty
#   cost = cost + sum(x %% 2 != 0) + penalty
#
#   # Add average and final state to history
#   history = c(history, list(avg), list(b))
#
#   # Return results (assuming the original function returns something)
#   return(list(
#     final_state = b,
#     cost = cost,
#     history = history,
#     ancestor = avg
#   ))
# }

find_greedy_distance_B_contig = function(a, b, penalty = 0) {
  # Make a copy of input vector to avoid modifying the original
  current = a

  # Calculate initial difference
  x = current - b

  # Initialize history and cost
  history = list(current)
  cost = 0

  # Continue until all differences are even (or zero)
  while (!all(x %% 2 == 0)) {
    #print(x)

    # Find positions with odd differences
    odd_positions = which(x %% 2 != 0)

    if (length(odd_positions) > 0) {
      # Find first position with odd difference
      head = min(odd_positions)

      # Find the end of the contiguous block (ignoring sign)
      end = head

      # Extend to consecutive positions with non-zero differences (regardless of sign)
      while (end < length(x) && x[end + 1] != 0) {
        end = end + 1
      }

      # Apply the operation based on the sign at the head position
      if (x[head] < 0) {
        current[head:end] = current[head:end] + 1
      } else {
        current[head:end] = current[head:end] - 1
      }

      # Update history and cost
      history = c(history, list(current))
      cost = cost + 1

      # Recalculate difference (keep all values, don't zero out even ones)
      x = current - b
    }
  }

  # Calculate average of final state and target
  r = rle(x != 0)
  avg = (current + b) / 2
  cost = cost + sum(r$values) + penalty

  avg = (current + b) / 2
  a_cost = cost - 1
  b_cost = 1

  # Add average and final state to history
  history = c(history, list(avg), list(b))

  return(list(
    history = history,
    cost = cost,
    ancestor = avg,
    b_cost = b_cost,
    a_cost = a_cost,
    a_delta = a - avg,
    b_delta = b - avg,
    left = current,
    right = b
  ))
}


# find_greedy_distance_with_avg_v1 = function(a, b, penalty = 0) {
#   # Make a copy of input vector to avoid modifying the original
#   current = a
#
#   # Calculate initial difference and set even differences to 0
#   x = current - b
#   cost = sum(rle(x %% 2 == 0)$value)
#   x[x %% 2 == 0] = 0
#
#   # Initialize history and cost
#   history = list(current)
#
#   # Continue until all differences are even (or zero)
#   while(!all(x %% 2 == 0)) {
#     # Find first position with non-zero difference
#     if(any(x != 0)) {
#       head = min(which(x != 0))
#
#       # Find the end of the contiguous block with same sign
#       end = head
#       sign_val = sign(x[head])
#
#       # Make sure not to go out of bounds
#       while(end < length(x) && sign(x[end + 1]) == sign_val) {
#         end = end + 1
#       }
#
#       # Apply the operation (add or subtract 1)
#       if(x[head] < 0) {
#         current[head:end] = current[head:end] + 1
#       } else {
#         current[head:end] = current[head:end] - 1
#       }
#
#       # Update history and cost
#       history = c(history, list(current))
#       cost = cost + 1
#
#       # Recalculate difference and set even differences to 0
#       x = current - b
#       x[x %% 2 == 0] = 0
#     }
#   }
#
#   # Calculate average of final state and target
#   avg = (current + b) / 2
#   cost = cost + 1 + penalty
#
#   # Add average and final state to history
#   history = c(history, list(avg), list(b))
#
#   return(list(
#     history = history,
#     cost = cost,
#     ancestor = avg
#   ))
# }
#
# find_greedy_distance_with_avg_v1 = function(a, b, penalty = 0) {
#   # Make a copy of input vector to avoid modifying the original
#   current = a
#
#   # Calculate initial difference and set even differences to 0
#   x = current - b
#   cost = sum(rle(x %% 2 == 0)$value)
#   x[x %% 2 == 0] = 0
#
#   # Initialize history and cost
#   history = list(current)
#
#   # Continue until all differences are even (or zero)
#   while(!all(x %% 2 == 0)) {
#     # Find first position with non-zero difference
#     if(any(x != 0)) {
#       head = min(which(x != 0))
#
#       # Find the end of the contiguous block with same sign
#       end = head
#       sign_val = sign(x[head])
#
#       # Make sure not to go out of bounds
#       while(end < length(x) && sign(x[end + 1]) == sign_val) {
#         end = end + 1
#       }
#
#       # Apply the operation (add or subtract 1)
#       if(x[head] < 0) {
#         current[head:end] = current[head:end] + 1
#       } else {
#         current[head:end] = current[head:end] - 1
#       }
#
#       # Update history and cost
#       history = c(history, list(current))
#       cost = cost + 1
#
#       # Recalculate difference and set even differences to 0
#       x = current - b
#       x[x %% 2 == 0] = 0
#     }
#   }
#
#   # Calculate average of final state and target
#   avg = (current + b) / 2
#   cost = cost + 1 + penalty
#
#   # Add average and final state to history
#   history = c(history, list(avg), list(b))
#
#   return(list(
#     history = history,
#     cost = cost,
#     ancestor = avg
#   ))
# }
#
# find_greedy_distance_with_avg_v2 = function(a, b, penalty = 0) {
#   # Make a copy of input vector to avoid modifying the original
#   current = a
#
#   # Calculate initial difference and set even differences to 0
#   x = current - b
#   cost = sum(x %% 2 == 0)
#   x[x %% 2 == 0] = 0
#
#   # Initialize history and cost
#   history = list(current)
#
#   # Continue until all differences are even (or zero)
#   while(!all(x %% 2 == 0)) {
#     # Find first position with non-zero difference
#     if(any(x != 0)) {
#       head = min(which(x != 0))
#
#       # Find the end of the contiguous block with same sign
#       end = head
#       sign_val = sign(x[head])
#
#       # Make sure not to go out of bounds
#       while(end < length(x) && sign(x[end + 1]) == sign_val) {
#         end = end + 1
#       }
#
#       # Apply the operation (add or subtract 1)
#       if(x[head] < 0) {
#         current[head:end] = current[head:end] + 1
#       } else {
#         current[head:end] = current[head:end] - 1
#       }
#
#       # Update history and cost
#       history = c(history, list(current))
#       cost = cost + 1
#
#       # Recalculate difference and set even differences to 0
#       x = current - b
#       x[x %% 2 == 0] = 0
#     }
#   }
#
#   # Calculate average of final state and target
#   avg = (current + b) / 2
#   cost = cost + 1 + penalty
#
#   # Add average and final state to history
#   history = c(history, list(avg), list(b))
#
#   return(list(
#     history = history,
#     cost = cost,
#     ancestor = avg
#   ))
# }
#
# find_greedy_distance_with_avg_v3 = function(a, b, penalty = 0) {
#   # Make a copy of input vector to avoid modifying the original
#   current = a
#
#   # Calculate initial difference and set even differences to 0
#   x = current - b
#   if (any(x %% 2 == 0)) {
#     cost = stats::median(abs(x[x %% 2 == 0]))
#     x[x %% 2 == 0] = 0
#   } else {
#     cost = 0
#   }
#
#   # Initialize history and cost
#   history = list(current)
#
#   # Continue until all differences are even (or zero)
#   while(!all(x %% 2 == 0)) {
#     # Find first position with non-zero difference
#     if(any(x != 0)) {
#       head = min(which(x != 0))
#
#       # Find the end of the contiguous block with same sign
#       end = head
#       sign_val = sign(x[head])
#
#       # Make sure not to go out of bounds
#       while(end < length(x) && sign(x[end + 1]) == sign_val) {
#         end = end + 1
#       }
#
#       # Apply the operation (add or subtract 1)
#       if(x[head] < 0) {
#         current[head:end] = current[head:end] + 1
#       } else {
#         current[head:end] = current[head:end] - 1
#       }
#
#       # Update history and cost
#       history = c(history, list(current))
#       cost = cost + 1
#
#       # Recalculate difference and set even differences to 0
#       x = current - b
#       x[x %% 2 == 0] = 0
#     }
#   }
#
#   # Calculate average of final state and target
#   avg = (current + b) / 2
#   cost = cost + 1 + penalty
#
#   # Add average and final state to history
#   history = c(history, list(avg), list(b))
#
#   return(list(
#     history = history,
#     cost = cost,
#     ancestor = avg
#   ))
# }
