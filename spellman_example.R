setwd("~/research/prelim")

# Hypothesis tests.

# Contained in base R.
COR <- "correlation"

# library(acepack)
# MAXCOR <- "maxcor"

library(energy)
DCOR <- "dcor"

library(dHSIC)
HSIC <- "hsic"

library(minerva)
# MIC <- "mic"
TIC <- "tic"

suppressPackageStartupMessages(library(HHG))
HHG <- "hhg"
NT = Fast.independence.test.nulltable(n = 23, nr.perm = 1000)

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

spellman <- read.csv(file = 'data/spellman_gene_expr_data.csv')
genes <- names(spellman)[2:ncol(spellman)]
x <- spellman[, "time"]

get_pvalues <- function(test_name) {
  apply_test <- function(gene) {
    y <- spellman[, gene]
    get_pvalue(x, y, test_name)
  }
  fname <- sprintf("results/spellman/%s_pvalues.txt", test_name)
  pvalues <- unlist(lapply(genes, apply_test))
  write.table(pvalues, file=fname, col.names=F, row.names=F)
}

# tests <- c(COR, DCOR, HSIC, HHG, TIC)
tests <- c(HHG)
for (test_name in tests) {
  get_pvalues(test_name)
}


