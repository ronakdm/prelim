# Contained in base R.
COR <- "correlation"
DCOR <- "dcor"
TIC <- "tic"
HHG <- "hhg"
set.seed(1)
NT = Fast.independence.test.nulltable(n = 23, nr.perm = 1000)


newcorrelation = function(x, y, pvalue = TRUE, discrete = TRUE, method = "asymptotic", iter = 1000, na.rm = FALSE, factor = FALSE) {
  # This is the main function that computes the new correlation coefficient and returns P-values. 
  if (na.rm == TRUE) {
    a = which(!is.na(x) & !is.na(y))
    x = x[a]
    y = y[a]
  }
  if (factor == TRUE) {
    if (!is.numeric(x)) x = as.numeric(factor(x))
    if (!is.numeric(y)) y = as.numeric(factor(y))
  }
  n = length(x)
  PI = rank(x, ties.method = "random")
  f = rank(y, ties.method = "max")/n
  g = rank(-y, ties.method = "max")/n
  ord = order(PI)
  f = f[ord]
  A1 = mean(abs(f[1:(n-1)] - f[2:n]))*(n-1)/(2*n)
  C = mean(g*(1-g))
  xi = 1 - A1/C
  if (pvalue == TRUE) {
    if (discrete == FALSE) return(list(xi = xi, pval = 1 - pnorm(sqrt(n)*xi/sqrt(2/5))))
    if (method == "asymptotic") {
      q = sort(f)
      ind = c(1:n)
      ind2 = 2*n - 2*ind + 1
      a = mean(ind2*q*q)/n
      c = mean(ind2*q)/n
      cq = cumsum(q)
      m = (cq + (n - ind)*q)/n
      b = mean(m^2)
      v = (a - 2*b + c^2)/(C^2)
      return(list(xi = xi, sd = sqrt(v), pval = 1 - pnorm(sqrt(n)*xi/sqrt(v))))
    }
    r = rep(0, iter)
    for (i in 1:iter) {
      x1 = runif(n, 0, 1)
      r[i] = newcorrelation(x1,y)$xi
    }
    return(list(xi = xi, sd = sqrt(var(r)*n), pval = mean(r >= xi)))
  }
  else return(xi)
}

fdr = function(pval, q = .05) {
  # Subroutine for calculating FDR.
  m = length(pval)
  l = c(1:m)*q/m
  ord = order(pval)
  u = pval[ord]
  w = which(u <= l)
  if (length(w) == 0) return(NULL)
  i = max(w)
  ind = ord[1:i]
  return(ind)
}

spellmantest = function(fdr = .05) {
  # Calculates P-values for the Spellman data for various tests of independence, and stores the P-values in a file.
  library(minerva)
  library(dHSIC)
  data(Spellman)
  library(HHG)
  library(energy)
  library(acepack)
  b = Spellman
  p = dim(Spellman)[2]
  n = dim(Spellman)[1]
  set.seed(1)
  NT = Fast.independence.test.nulltable(n = n, nr.perm = 1000) 
  val.cor = val.xi = val.dcor = val.mic = val.maxcor = rep(0, p-1)
  p.cor = p.xi = p.dcor = p.tic = p.hhg = p.hsic = rep(0, p-1)
  for (i in 2:p) {
    print(i)
    val.cor[i-1] = cor(b[,1], b[,i])
    p.cor[i-1] = cor.test(b[,1], b[,i])$p.value
    a = newcorrelation(b[,1], b[,i], discrete = TRUE)
    val.xi[i-1] = a$xi
    p.xi[i-1] = a$pval
    val.dcor[i-1] = dcor(b[,1], b[,i])
    p.dcor[i-1] = dcor.test(b[,1], b[,i], R = 1000)$p.value
    val.mic[i-1] = cstats(as.matrix(b[,1]), as.matrix(b[,i]))[3]
    a = ace(b[,1], b[,i])
    val.maxcor[i-1] = cor(a$tx, a$ty)
    p.hhg[i-1] = Fast.independence.test(b[,1], b[,i], NullTable = NT, combining.type = 'Fisher')$Fisher.pvalue
    p.tic[i-1] = as.numeric(mictools(cbind(b[,1], b[,i]), nperm = 1000)$pval[1])
    p.hsic[i-1] = dhsic.test(b[,1], b[,i])$p.value
  }
  ind.cor = fdr(p.cor, q = fdr)+1
  ind.xi = fdr(p.xi, q = fdr)+1
  ind.dcor = fdr(p.dcor, q = fdr)+1
  ind.tic = fdr(p.tic, q = fdr)+1
  ind.hhg = fdr(p.hhg, q = fdr)+1
  ind.hsic = fdr(p.hsic, q = fdr)+1
  s = list(val.cor = val.cor, val.dcor = val.dcor, val.xi =val.xi, val.mic = val.mic, val.maxcor = val.maxcor, p.cor = p.cor, p.xi = p.xi, p.dcor = p.dcor, p.tic = p.tic, p.hhg = p.hhg, p.hsic = p.hsic, ind.cor = ind.cor, ind.dcor = ind.dcor, ind.xi = ind.xi, ind.tic = ind.tic, ind.hhg = ind.hhg, ind.hsic = ind.hsic)
  save(s, file = "/Users/sourav/Desktop/Tex Files/correlation-testing/spellmantest.RData")
  return(s)
}

spellmantest_individual = function(fdr = .05) {
  # Calculates P-values for the Spellman data for various tests of independence, and stores the P-values in a file.
  library(minerva)
  data(Spellman)
  b = Spellman
  p = dim(Spellman)[2]
  n = dim(Spellman)[1]
  stat_vals = rep(0, p-1)
  p_vals = rep(0, p-1)
  library(HHG)
  for (i in 2:p) {
    print(i)
    # stat_vals[i-1] = dcor(b[,1], b[,i])
    set.seed(1)
    p_vals[i-1] = Fast.independence.test(b[,1], b[,i], NullTable = NT, combining.type = 'Fisher')$Fisher.pvalue
  }
  inds = fdr(p_vals, q = fdr)
  print(length(inds))
  
  test_name <- HHG
  prefix <- sprintf("results/spellman/%s_reference_", test_name)
  write.table(unlist(stat_vals), file=paste(prefix, "stats.txt", sep=""), col.names=F, row.names=F)
  write.table(unlist(p_vals), file=paste(prefix, "pvals.txt", sep=""), col.names=F, row.names=F)
  write.table(unlist(inds), file=paste(prefix, "inds.txt", sep=""), col.names=F, row.names=F)
}

spellmantest_selection = function(fdr = .05) {
  # Calculates P-values for the Spellman data for various tests of independence, and stores the P-values in a file.
  data(Spellman)
  b = Spellman
  p = dim(Spellman)[2]
  n = dim(Spellman)[1]
  set.seed(1)
  NT = Fast.independence.test.nulltable(n = n, nr.perm = 1000) 
  p.cor = p.xi = p.dcor = p.tic = p.hhg = p.hsic = rep(0, p-1)
  # for (i in 2:10) {
  for (i in 2:p) {
    print(i)
    # p.cor[i-1] = cor.test(b[,1], b[,i])$p.value
    # a = newcorrelation(b[,1], b[,i], discrete = TRUE)
    # p.xi[i-1] = a$pval
    # set.seed(1)
    # p.dcor[i-1] = dcor.test(b[,1], b[,i], R = 1000)$p.value
    # set.seed(1)
    p.hhg[i-1] = Fast.independence.test(b[,1], b[,i], NullTable = NT, combining.type = 'Fisher')$Fisher.pvalue
    # set.seed(1)
    # p.tic[i-1] = as.numeric(mictools(cbind(b[,1], b[,i]), nperm = 1000)$pval[1])
  }
  # ind.cor = fdr(p.cor, q = fdr)+1
  # print(length(ind.cor))
  # ind.xi = fdr(p.xi, q = fdr)+1
  # print(length(ind.xi))
  # ind.dcor = fdr(p.dcor, q = fdr)+1
  # print(length(ind.dcor))
  # ind.tic = fdr(p.tic, q = fdr)+1
  # print(length(ind.tic))
  ind.hhg = fdr(p.hhg, q = fdr)+1
  print(length(ind.hhg))
  
  # xionly = setdiff(xi, union(cor, union(dcor, union(tic, union(hhg, hsic)))))
  # notxi = setdiff(union(cor, union(dcor, union(tic, union(hhg, hsic)))), xi)
  # xionly = setdiff(ind.xi, union(ind.cor, union(ind.dcor, union(ind.tic, ind.hhg))))
  # notxi = setdiff(union(ind.cor, union(ind.dcor, union(ind.tic, ind.hhg))), ind.xi)
  # prefix <- "results/spellman/full_reference_"
  # write.table(unlist(xionly), file=paste(prefix, "xionly.txt", sep=""), col.names=F, row.names=F)
  # write.table(unlist(notxi), file=paste(prefix, "notxi.txt", sep=""), col.names=F, row.names=F)
}

# library(minerva)
# library(dHSIC)
# library(HHG)
# library(energy)
# library(acepack)

spellmantest_selection()