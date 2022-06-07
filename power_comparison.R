setwd("~/research/prelim")
COR <- "correlation"

# Hyperparameters.
NOISE_LEVELS <- 10
N <- 100

# Hypothesis tests.
library(energy)
DCOR <- "dcor"

library(dHSIC)
HSIC <- "hsic"

library(minerva)
MIC <- "mic"
TIC <- "tic"

suppressPackageStartupMessages(library(HHG))
HHG <- "hhg"
NT <- Fast.independence.test.nulltable(n = N)

# Relationships.
LINEAR <- "linear"
STEP_FUNC <- "step_function"
W_SHAPED <- "w_shaped"
SINUSOID <- "sinusoid"
CIRCULAR <- "circular"
HETERO <- "heteroskedastic"

get_pvalue <- function(x, y, test_name) {
  if (test_name == DCOR) {
    set.seed(1)
    return(dcor.test(x, y, R = 1000)$p.value)
  } else if (test_name == HSIC) {
    set.seed(1)
    return(dhsic.test(x, y)$p.value)
  } else if (test_name == HHG) {
    set.seed(1)
    return(Fast.independence.test(x, y, NullTable = NT, combining.type = "Fisher")$Fisher.pvalue)
  } else if (test_name == TIC){
    set.seed(1)
    return(as.numeric(mictools(cbind(x, y), nperm = 1000)$pval[1]))
  } else if (test_name == COR) {
    return(cor.test(x, y)$p.value)
  } else {
    stop("Unrecognized 'test_name': ", test_name)
  }
}


run_test <- function(x, y, test_name) {
  if (test_name == DCOR) {
    # from energy package.
    return(dcor(x, y))
  } else if (test_name == HSIC) {
    return(dhsic(x, y)$dHSIC)
  } else if (test_name == MIC) {
    return(cstats(as.matrix(x), as.matrix(y))[3])
  } else if (test_name == HHG) {
    return(Fast.independence.test(x, y, NullTable = NT, combining.type = "Fisher")$Fisher)
  } else {
    stop("Unrecognized 'test_name': ", test_name)
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

get_true_power <- function(test_name, relationship, noise_level) {
  # critical_value <- get_critical_value(test_name, relationship, noise_level)
  fname <- sprintf("data/power/%s_noise_level_%d_joint.csv", relationship, noise_level)
  df <- read.csv(file = fname)
  n_col <- ncol(df)
  apply_test <- function(i) {
    x <- df[, 2 * i - 1]
    y <- df[, 2 * i]
    pval <- get_pvalue(x, y, test_name)
    return(as.integer(pval <= 0.05))
  }
  rejects <- unlist(lapply(1:(ncol(df) / 2), apply_test))
  power <- mean(rejects)
  power
}

get_powers <- function(test_name, relationship) {
  powers <- c()
  for (noise_level in 0:NOISE_LEVELS) {
    powers <- c(powers, get_power(test_name, relationship, noise_level))
  }
  fname <- sprintf("results/power/%s_power_%s.txt", test_name, relationship)
  write.table(unlist(powers), file = fname, col.names = F, row.names = F)
}

get_true_powers <- function(test_name, relationship) {
  powers <- c()
  for (noise_level in 0:NOISE_LEVELS) {
    powers <- c(powers, get_true_power(test_name, relationship, noise_level))
  }
  fname <- sprintf("results/power/%s_true_power_%s.txt", test_name, relationship)
  write.table(unlist(powers), file = fname, col.names = F, row.names = F)
}



relationships <- c(STEP_FUNC, W_SHAPED, SINUSOID, CIRCULAR, HETERO)
# relationships <- c(LINEAR)
# test_names <- c(DCOR, HSIC, MIC, HHG)
test_names <- c(COR)
for (test_name in test_names) {
  for (relationship in relationships) {
    print(sprintf("Computing '%s' power on '%s' relationship...", test_name, relationship))
    get_true_powers(test_name, relationship)
  }
}
