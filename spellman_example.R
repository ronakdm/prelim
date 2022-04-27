setwd("~/research/prelim")

# Data.
N <- 100

# Hypothesis tests.

# Contained in base R.
COR <- "correlation"

library(acepack)
MAXCOR <- "maxcor"

library(energy)
DCOR <- "dcor"

library(dHSIC)
HSIC <- "hsic"

library(minerva)
MIC <- "mic"
TIC <- "tic"

library(HHG)
HHG <- "hhg"
NT = Fast.independence.test.nulltable(n = N)

# Get p-value of corresponding test. No MAXCOR and MIC from paper.
get_pvalue <- function(x, y, test_name) {
  if (test_name == DCOR) {
    return(dcor.test(x, y, R = 1000)$p.value)
  } else if (test_name == HSIC) {
    return(dhsic.test(x, y)$p.value)
  } else if (test_name == HHG) {
    return(Fast.independence.test(x,y, NullTable = NT, combining.type = 'Fisher')$Fisher.pvalue)
  } else if (test_name == TIC){
    return(as.numeric(mictools(cbind(x, y), nperm = 1000)$pval[1]))
  } else if (test_name == COR) {
    return(cor.test(x, y)$p.value)
  } else {
    throw("Unrecognized 'test_name': ", test_name)
  }
}