library(energy)

# Hyperparameters.
NOISE_LEVELS <- 10

# Relationships.
LINEAR <- "linear"
STEP_FUNC <- "step_function"
W_SHAPED <- "w_shaped"
SINUSOID <- "sinusoid"
CIRCULAR <- "circular"
HETERO <- "heteroskedastic"

DCOR <- "dcor"

run_test <- function(x, y, test_name) {
  if (test_name == DCOR) {
    # from energy package.
    return(dcor(x, y))
  }
}

get_critical_value <- function(test_name, relationship, noise_level) {
  fname <- sprintf("data/power/%s_noise_level_%d_marginal.csv", relationship, noise_level)
  df <- read.csv(file = fname)
  n_col <- ncol(df)
  apply_test <- function(i) {
    x <- df[, 2 * i - 1]
    y <- df[, 2 * i]
    run_test(x, y, test_name)
  }
  test_stat_dist <- unlist(lapply(1:(ncol(df) / 2), apply_test))
  critical_value <- quantile(test_stat_dist, 0.95)
  critical_value
}

get_power <- function(test_name, relationship, noise_level) {
  critical_value <- get_critical_value(test_name, relationship, noise_level)
  fname <- sprintf("data/power/%s_noise_level_%d_joint.csv", relationship, noise_level)
  df <- read.csv(file = fname)
  n_col <- ncol(df)
  apply_test <- function(i) {
    x <- df[, 2 * i - 1]
    y <- df[, 2 * i]
    run_test(x, y, test_name)
  }
  rejects <- (unlist(lapply(1:(ncol(df) / 2), apply_test)) > critical_value)
  power <- mean(rejects)
  power
}

get_powers <- function(test_name, relationship) {
  powers <- c()
  for (noise_level in 0:NOISE_LEVELS ) {
    powers <- c(powers, get_power(test_name, relationship, noise_level))
  }
  fname <- sprintf("results/power/%s_power_%s.txt", test_name, relationship)
  write.table(unlist(powers), file=fname, col.names=F, row.names=F)
}

test_name <- DCOR
relationship <- CIRCULAR
noise_level <- 0
