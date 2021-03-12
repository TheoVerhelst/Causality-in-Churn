# Derives the parameters of the bivariate beta distribution such that
# the resulting population respects the proportions beta, gamma and delta.
derive_bivariate_beta_parameters <- function(beta, gamma, delta, A = 1) {
  alpha <- 1 - beta - gamma - delta
  v <- beta * gamma - delta * alpha
  a_00 <- -v + A * alpha
  a_01 <-  v + A * gamma
  a_10 <-  v + A * beta
  a_11 <- -v + A * delta

  cov_S0_S1 <- (a_11*a_00 - a_10*a_01) / (A * (A + 1))
  E_S0 <- (a_11 + a_10) / A
  E_S1 <- (a_11 + a_01) / A

  return(list(
    a = c(a_00, a_01, a_10, a_11),
    cov_S0_S1 = cov_S0_S1,
    E_S0 = E_S0,
    E_S1 = E_S1
  ))
}

# Create N samples from a Dirichlet distribution with parameter vector alpha.
sample_dirichlet <- function(N, alpha) {
  stopifnot(all(alpha > 0))
  gamma_sample <- sapply(alpha, function(alpha_i)
    rgamma(N, shape = alpha_i, rate = 1)
  )
  # Transpose the result so we have N rows and length(alpha) columns
  return(aperm(sapply(1:N, function(j) gamma_sample[j,] / sum(gamma_sample[j,]))))
}

# Creates N samples from a bivariate beta distribution described in
# "Olkin, Ingram, and Thomas A. Trikalinos. "Constructions for a bivariate
# beta distribution." Statistics & Probability Letters 96 (2015): 54-60."
# The treatment and target variables are then sampled from Bernoulli
# distributions. Beta, gamma and delta dictate the distribution of Y_0 and Y_1.
# noise_S_0 and noise_S_1 are (std dev of) Gaussian noise added to S_0 and S_1.
# A controls the spread of the beta marginals.
# use_churn_convention determines which value of Y_0, Y_1 is considered
# persuadable.
sample_bivariate_beta <- function(beta, gamma, delta,
                                  N = 20000, proba_treatment = 0.65, A = 1,
                                  noise_S_0 = 0.01, noise_S_1 = 0.01,
                                  use_churn_convention = TRUE) {
  parameters <- derive_bivariate_beta_parameters(beta, gamma, delta, A = A)
  U <- sample_dirichlet(N=N, alpha=parameters$a)

  data <- data.frame(
    S_0 = U[,4] + U[,3],
    S_1 = U[,4] + U[,2]
  )

  if (use_churn_convention) {
    data$uplift <- data$S_0 - data$S_1
  } else {
    data$uplift <- data$S_1 - data$S_0
  }


  data$Y_0 <- as.factor(as.numeric(runif(N) <= data$S_0))
  data$Y_1 <- as.factor(as.numeric(runif(N) <= data$S_1))

  data$T_ <- as.factor(as.numeric(runif(N) <= proba_treatment))
  data$Y <- as.factor(ifelse(
    data$T_ == "0",
    as.integer(as.character(data$Y_0)),
    as.integer(as.character(data$Y_1))
  ))

  # Add random noise (increase variance) of the scores
  data$S_0_hat <- data$S_0 + rnorm(n = N, sd = noise_S_0)
  data$S_1_hat <- data$S_1 + rnorm(n = N, sd = noise_S_1)

  # Compute the estimated uplift
  if (use_churn_convention) {
      data$uplift_hat <- data$S_0_hat - data$S_1_hat
    } else {
      data$uplift_hat <- data$S_1_hat - data$S_0_hat
  }

  return(data)
}


# Creates N samples from a unit simplex of dimension k
# Source: Smith, Noah A., and Roy W. Tromble. "Sampling uniformly from
# the unit simplex." Johns Hopkins University, Tech. Rep 29 (2004).
sample_simplex <- function(N, k) {
  aperm(sapply(1:N, function(i) {
    # Sample k - 1 reals from [0, 1] without replacement
    numbers <- unique(runif(n = k - 1))
    remaining <- k - 1 - length(numbers)
    while (remaining > 0) {
      numbers <- unique(c(numbers, runif(n = remaining)))
      remaining <- k - 1 - length(numbers)
    }

    numbers <- c(0, sort(numbers), 1)
    return(sapply(1:k, function(i) numbers[i + 1] - numbers[i]))
  }))
}
