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


permtest_HHG = function(n, iter = 200) {
# This function returns the computational time for the HHG test of independence in a sample of size n, based on a permutation test of size iter.
	library(HHG)
	x = rnorm(n, 0, 1)
	y = rnorm(n, 0, 1)
	s0 = Sys.time()
	Dy = as.matrix(dist(y))
	Dx = as.matrix(dist(x))
	h = hhg.test(Dx, Dy, nr.perm = 200)
	s1 = Sys.time()
	print(s1-s0)	
}


permtest_MIC = function(n, iter = 200) {
# This function returns the computational time for the MIC test of independence in a sample of size n, based on a permutation test of size iter.
	library(minerva)
	x = rnorm(n, 0, 1)
	y = rnorm(n, 0, 1)
	s0 = Sys.time()
	for (i in 1:iter) {
		x1 = sample(x, size = n, replace = F)
		m = cstats(as.matrix(x1),as.matrix(y))
	}
	s1 = Sys.time()
	print(s1-s0)
}

permtest_HSIC = function(n, iter = 200) {
# This function returns the computational time for the HSIC test of independence in a sample of size n, based on a permutation test of size iter.
	library(dHSIC)
	x = rnorm(n, 0, 1)
	y = rnorm(n, 0, 1)
	s = system.time(a <- dhsic.test(x,y, B = iter))
	print(s)
}


permtest_dcor = function(n, iter = 200) {
# This function returns the computational time for the dCor test of independence in a sample of size n, based on a permutation test of size iter.
	library(energy)
	x = rnorm(n, 0, 1)
	y = rnorm(n, 0, 1)
	s = system.time(a <- dcor.test(x,y, R = iter))
	print(s)
}

permtest_newcorrelation = function(n, iter = 200) {
# This function returns the computational time for the test of independence based on the new correlation coefficient in a sample of size n, using the asymptotic theory for the test.
	x = rnorm(n, 0, 1)
	y = rnorm(n, 0, 1)
	s = system.time(a <- newcorrelation(x,y, pval = TRUE))
	print(s)
}



xihist = function(n, iter, bins = 50) {
# Draws the histogram of xi values versus theoretical null distribution for continuous random variables. 
	r = rep(0, iter)
	for (i in 1:iter) {
		x = runif(n, 0, 1)
		y = runif(n, 0, 1)
		r[i] = newcorrelation(x,y)$xi
	}
	r = r*sqrt(n)
	hist(r, freq = F, breaks = bins, main = "", xlab = "", ylab = "", xlim = c(-3.5,3.5), ylim = c(0, .75), col = "light grey")
	x = seq(min(r), max(r), length = 500)
	y = dnorm(x, 0, sqrt(2/5))
	lines(x,y)
}

xihist2 = function(n, iter, bins = 50) {
# Draws the histogram of xi values versus theoretical null distribution for discrete random variables. 
	r = rep(0, iter)
	s = rep(0, iter)
	for (i in 1:iter) {
		x = rbinom(n, 3, .5)
		y = rbinom(n, 3, .5)
		a = newcorrelation(x, y, discrete = T)
		r[i] = a$xi
		s[i] = a$sd
	}
	sd = a$sd
	r = r*sqrt(n)
	hist(r, freq = F, breaks = bins, main = "", xlab = "", ylab = "", xlim = c(-3.5,3.5), ylim = c(0, .75), col = "light grey")
	x = seq(min(r), max(r), length = 500)
	y = dnorm(x, 0, sd = sd)
	lines(x,y)
	print(c(min(s), max(s)))
	print(c(quantile(s,.025), quantile(s,.975)))
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

spellmangraph = function(s) {
# Using the P-values for tests of independence in the Spellman data obtained using the spellmantest function, sorts the variables by order of importance. Produces the graphs displayed in the paper and stores them in separate pdf files.
	library(minerva)
	library(caret)
	data(Spellman)
	b = Spellman
	names = colnames(Spellman)
	xi = s$ind.xi
	cor = s$ind.cor 
	dcor = s$ind.dcor
	tic = s$ind.tic
	hhg = s$ind.hhg
	hsic = s$ind.hhg
	pxi = s$p.xi
	pcor = s$p.cor
	pdcor = s$p.dcor
	ptic = s$p.tic
	phhg = s$p.hhg 
	phsic = s$p.hhg
	xionly = setdiff(xi, union(cor, union(dcor, union(tic, union(hhg, hsic)))))
	notxi = setdiff(union(cor, union(dcor, union(tic, union(hhg, hsic)))), xi)
	notxicor = pcor[notxi-1]
	notxidcor = pdcor[notxi-1]
	notxitic = ptic[notxi-1]
	notxihhg = phhg[notxi-1]
	notxihsic = phsic[notxi-1]
	p = pmin(notxicor, notxidcor, notxitic, notxihhg, notxihsic)
	u = order(p)
	notxiord = notxi[u]
	ind1 = xionly[1:6]
	ind2 = notxiord[1:6]
	ind3 = sample(xionly, size = 6, replace = F)
	ind4 = sample(notxi, size = 6, replace = F)
	pdf("/Users/sourav/Desktop/Tex Files/correlation-testing/newcorrelation-spellman1.pdf")
# This graph displays the top six variables selected by xi.
	par(mfrow = c(3,2), cex = .5)
	for (i in ind1) {
		x = b[,1]
		y = b[,i]
		plot(x,y, ylab="", xlab=names[i], xaxt = 'n', yaxt = 'n')
		fit = knnreg(as.matrix(x), y, k=3)
		lines(x, predict(fit, x), lty = "dashed")
	}
	dev.off()
	pdf("/Users/sourav/Desktop/Tex Files/correlation-testing/newcorrelation-spellman2.pdf")
# This graph displays the top six variables selected by other coefficients.
	par(mfrow = c(3,2), cex = .5)
	for (i in ind2) {
		x = b[,1]
		y = b[,i]
		plot(x,y, ylab="", xlab=names[i], xaxt = 'n', yaxt = 'n')
		fit = knnreg(as.matrix(x), y, k=3)
		lines(x, predict(fit, x), lty = "dashed")
	}
	dev.off()
	pdf("/Users/sourav/Desktop/Tex Files/correlation-testing/newcorrelation-spellman3.pdf")
# This graph displays six randomly chosen variables among those selected by xi.
	par(mfrow = c(3,2), cex = .5)
	for (i in ind3) {
		x = b[,1]
		y = b[,i]
		plot(x,y, ylab="", xlab=names[i], xaxt = 'n', yaxt = 'n')
		fit = knnreg(as.matrix(x), y, k=3)
		lines(x, predict(fit, x), lty = "dashed")
	}
	dev.off()
	pdf("/Users/sourav/Desktop/Tex Files/correlation-testing/newcorrelation-spellman4.pdf")
# This graph displays six randomly chosen variables among those selected by other variables.
	par(mfrow = c(3,2), cex = .5)
	for (i in ind4) {
		x = b[,1]
		y = b[,i]
		plot(x,y, ylab="", xlab=names[i], xaxt = 'n', yaxt = 'n')
		fit = knnreg(as.matrix(x), y, k=3)
		lines(x, predict(fit, x), lty = "dashed")
	}
	dev.off()
}

hist_dep = function(n, iter, bins = 50) {
# Draws the histogram of xi values for dependent data.
	r = rep(0, iter)
	for (i in 1:iter) {
		x = rbinom(n, 1, .4)
		z = rbinom(n, 1, .5)
		y = x*z
		r[i] = newcorrelation(x,y)$xi
	}
	hist(r, freq = F, breaks = bins, main = "", xlab = "", ylab = "",  col = "light grey")
	m = mean(r)
	s = sd(r)
	x = seq(min(r), max(r), length = 500)
	y = dnorm(x, m, s)
	lines(x,y)
}
