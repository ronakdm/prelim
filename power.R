library(energy)

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
    return(dcor(x, y))
  }
}

n <- 100
x <- rnorm(n)
y <- rnorm(n)
result <- run_test(x, y, DCOR)

