uplift_curve <- function(Y, T_, score, use_churn_convention = TRUE) {
  stopifnot(length(Y) == length(T_))
  stopifnot(length(T_) == length(score))

  N <- length(Y)

  score_ranking <- order(score, decreasing = TRUE)
  Y <- as.numeric(as.character(Y[score_ranking]))
  T_ <- as.numeric(as.character(T_[score_ranking]))

  grid <- data.frame(n = 1:N,
                     normalized_x = (1:N) / N,
                     n_1 = cumsum(T_),
                     n_0 = cumsum(1 - T_),
                     response_1 = cumsum(Y & T_),
                     response_0 = cumsum(Y & (1 - T_)))
  if(use_churn_convention) {
    grid$normalized_uplift <- (grid$response_0 / grid$n_0 - grid$response_1 / grid$n_1) * grid$n / N
  } else {
    grid$normalized_uplift <- (grid$response_1 / grid$n_1 - grid$response_0 / grid$n_0) * grid$n / N
  }

  return(grid)
}

# Computes the uplift curve assuming we know the true underlying uplift for each
# customers. This is different from the empirical uplift curve above where
# the uplift is locally estimated using target and control groups.
true_uplift_curve <- function(score, uplift, use_churn_convention = TRUE) {
  stopifnot(length(score) == length(uplift))

  N <- length(score)

  uplift <- uplift[order(score, decreasing = TRUE)]
  grid <- data.frame(n = 1:N,
                     normalized_x = (1:N) / N,
                     normalized_uplift = cumsum(uplift) * (1:N) / N)

  return(grid)
}


AUUC <- function(Y, T_, score, use_churn_convention = TRUE) {
  grid <- uplift_curve(Y, T_, score, use_churn_convention = use_churn_convention)
  return(mean(grid$normalized_uplift, na.rm = TRUE))
}

true_AUUC <- function(score, uplift, use_churn_convention = TRUE) {
  grid <- true_uplift_curve(score, uplift, use_churn_convention = use_churn_convention)
  return(mean(grid$normalized_uplift, na.rm = TRUE))
}

# Computes the causal precision curve, by supposing that we
# know the outcome conditional probabilities S_0 = P(Y_0 = 1 | X) and S_1 = P(Y_1 = 1 | X).
# We don't even need the value of Y_0 or Y_1, since these are determined by S_0 and S_1.
causal_precision_curve <- function(score, S_0, S_1, n_steps = 100, use_churn_convention = TRUE) {
  stopifnot(length(score) == length(S_0))
  stopifnot(length(S_0) == length(S_1))

  threshold_axis <- (5:n_steps) / n_steps
  ranking <- order(score, decreasing = TRUE)

  return(do.call(rbind, lapply(threshold_axis, function(threshold) {
    n_included_custs <- as.integer(threshold * length(score))
    tau <- score[ranking[n_included_custs]]
    predicted_persuadable <- score > tau

    if (use_churn_convention) {
      precision <- mean(S_0[predicted_persuadable] * (1 - S_1[predicted_persuadable]))
    } else {
      precision <- mean((1 - S_0[predicted_persuadable]) * S_1[predicted_persuadable])
    }

    return(data.frame(tau = tau, precision = precision, threshold = threshold))
  })))
}


AUCPC <- function(score, S_0, S_1, n_steps = 100, use_churn_convention = TRUE) {
  return(mean(causal_precision_curve(score, S_0, S_1, n_steps, use_churn_convention)$precision, na.rm = TRUE))
}


# Calibration of posterior probabilities as shown in
# Dal Pozzolo, Andrea, et al. "Calibrating probability
# with undersampling for unbalanced classification."
# 2015 IEEE Symposium Series on Computational Intelligence.
# IEEE, 2015.
#
# The scores and Y do not have to have the same length.
# In particular, we can calibrate scores for which we don't
# know Y (such as calibrating S_0 for customers where T = 1).
calibrate_score <- function(score, Y) {
  P_Y_1 <- mean(Y == "1")
  P_Y_0 <- mean(Y == "0")
  ratio <- P_Y_1 / P_Y_0
  return(ratio * score / (ratio * score - score + 1))
}
