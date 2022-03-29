library(energy)
library(minerva)
library(HHG)
library(acepack)
library(dHSIC)

nsim=500                          # The number of null datasets we use to estimate our rejection reject regions for an alternative with level 0.05
nsim2 = 500                           # Number of alternative datasets we use to estimate our power

num.noise <- 10                  # The number of different noise levels used
noise <- 3                        # A constant to determine the amount of noise

n=100                          # Number of data points per simulation

val.hsic = val.new=val.cor=val.dcor=val.mic=val.tic=val.maxcor=val.hhg =rep(NA,nsim)              # Vectors holding the null "correlations" (for the various correlation coefficients that we are comparing) for each of the nsim null datasets at a given noise level

val.hsic2 = val.new2=val.cor2=val.dcor2=val.mic2=val.tic2=val.maxcor2=val.hhg2= rep(NA,nsim2)              # Vectors holding the alternative "correlations" (for the various correlation coefficients that we are comparing) for each of the nsim2 alternative datasets at a given noise level

power.hsic = power.new=power.cor=power.dcor=power.mic=power.tic= power.maxcor=power.hhg = array(NA, c(6,num.noise+1))                # Arrays holding the estimated power for each of the "correlation" types, for each data type (linear, parabolic, etc...) with each noise level

## We loop through the noise level and functional form; each time we estimate a null distribution based on the marginals of the data, and then use that null distribution to estimate power

## We use a uniformly distributed x.

NT = Fast.independence.test.nulltable(n = n)


for(l in 0:num.noise){
	print(l)
	ll = l+1
  for(typ in 1:6){

    ## This next loop simulates data under the null with the correct marginals (x is uniform, and y is a function of a uniform with gaussian noise)
    
    for(ii in 1:nsim){
      x=runif(n,-1,1)
      								#Linear + noise
      if(typ==1){
        y=.5*x+ noise *(l/num.noise)* rnorm(n)
      }
									# Step function + noise
      if(typ==2){
		 y = -3*(x< -.5) + 2*(x >= -.5 & x<0) - 4*(x>=0 & x < .5) - 3*(x>=.5) + 10*(l/num.noise)*rnorm(n)
      }

									#W-shape + noise
 	  if (typ == 3) {
  	  	y = abs(x+.5)*(x<0) + abs(x-.5)*(x>=0) + .25*noise * (l/num.noise) *rnorm(n)
 	  }

					
                                              #cosine + noise
      if(typ==4){
        y=cos(8*pi*x) + noise * (l/num.noise) *rnorm(n)
      }


        
      if (typ == 5) {               #circle + noise
      	y = (2*rbinom(n,1,.5)-1)*sqrt(1 - x^2) + .3*noise * (l/num.noise) *rnorm(n)
      }


      if (typ == 6) {               #heteroskedastic
      	y =  noise *  (3*(abs(x)<.5)*(1-l/num.noise) + l/num.noise) *rnorm(n)
      }

      
      x <- runif(n,-1,1)                       # We resimulate x so that we have the null scenario
      
      val.cor[ii]=(cor(x,y))^2           
      val.dcor[ii]=dcor(x,y)             
      u = cstats(as.matrix(x),as.matrix(y))
      a = u[3]
      b = u[4]
      val.mic[ii]= a          
      val.hhg[ii] = Fast.independence.test(x,y, NullTable = NT, combining.type = 'Fisher')$Fisher
	  val.new[ii]=newcorrelation(x,y)$xi
	  val.hsic[ii] = dhsic(x,y)$dHSIC
    }


    ## Next we calculate our rejection cutoffs
    
    cut.cor=quantile(val.cor,.95)
    cut.dcor=quantile(val.dcor,.95)
    cut.mic=quantile(val.mic,.95)
    cut.hhg=quantile(val.hhg,.95)
    cut.new=quantile(val.new,.95)
    cut.hsic = quantile(val.hsic, .95)
    

    ## Next we simulate the data again, this time under the alternative
    
    for(ii in 1:nsim2){
      x=runif(n,-1,1)

      								#Linear + noise
      if(typ==1){
        y=.5*x+ noise *(l/num.noise)* rnorm(n)
      }

									# Step function + noise
      if(typ==2){
		 y = -3*(x< -.5) + 2*(x >= -.5 & x<0) - 4*(x>=0 & x < .5) -1*(x>=.5) + 10*(l/num.noise)*rnorm(n)
     }

									#W-shape + noise
 	  if (typ == 3) {
  	  	y = abs(x+.5)*(x<0) + abs(x-.5)*(x>=0) + .25*noise * (l/num.noise) *rnorm(n)
 	  }

                                              # Cosine + noise
      if(typ==4){
        y=cos(8*pi*x) + noise * (l/num.noise) *rnorm(n)
      }

      if (typ == 5) {               # circle + noise
      	y = (2*rbinom(n,1,.5)-1)*sqrt(1 - x^2) + .3*noise * (l/num.noise) *rnorm(n)
      }


      if (typ == 6) {               #heteroskedastic
      	y =  noise * (3*(abs(x)<.5)*(1-l/num.noise) + l/num.noise)*rnorm(n)
      }


      ## We again calculate our "correlations"
      
      val.cor2[ii]=(cor(x,y))^2
      val.dcor2[ii]=dcor(x,y)
      u = cstats(as.matrix(x),as.matrix(y))
      a = u[3]
      b = u[4]
      val.mic2[ii]= a          
      val.hhg2[ii] = Fast.independence.test(x,y, NullTable=NT, combining.type = 'Fisher')$Fisher
	  val.new2[ii]=newcorrelation(x,y)$xi
 	  val.hsic2[ii] = dhsic(x,y)$dHSIC     
    }

    ## Now we estimate the power as the number of alternative statistics exceeding our estimated cutoffs
    
    power.cor[typ,ll] <- sum(val.cor2 > cut.cor)/nsim2
    power.dcor[typ,ll] <- sum(val.dcor2 > cut.dcor)/nsim2
    power.mic[typ,ll] <- sum(val.mic2 > cut.mic)/nsim2
    power.hhg[typ,ll] <- sum(val.hhg2 > cut.hhg)/nsim2
    power.new[typ,ll] <- sum(val.new2 > cut.new)/nsim2
    power.hsic[typ,ll] <- sum(val.hsic2 > cut.hsic)/nsim2
    
  }
}


noise.levels = (0:num.noise)/num.noise

save.image()

## The rest of the code is for plotting the image


pdf("/Users/sourav/Desktop/Tex Files/correlation-testing/newcorrelation-powerplot9.pdf")
#par(mfrow = c(4,2), cex = 0.45)
# par(mfrow = c(2,2), cex = .7)
par(mfrow = c(3,2), cex = 0.6)
plot(noise.levels, power.dcor[1,], ylim = c(0,1), xlim = c(0,1), main = "Linear", xlab = "Noise Level", ylab = "Power", pch = 1, type = 'b')
points(noise.levels, power.mic[1,], pch = 2, type = 'b')
points(noise.levels, power.hsic[1,], pch = 3, type = 'b')
points(noise.levels, power.hhg[1,], pch = 4, type = 'b')
points(noise.levels, power.new[1,], pch = 5, type = 'b')
legend("topright",c("dcor", "MIC", "HSIC", "HHG", expression(xi[n])), pch = c(1,2,3,4,5))

plot(noise.levels, power.dcor[2,], ylim = c(0,1), xlim = c(0,1),  main = "Step Function", xlab = "Noise Level", ylab = "Power", pch = 1, type = 'b')
points(noise.levels, power.mic[2,], pch = 2, type = 'b')
points(noise.levels, power.hsic[2,], pch = 3, type = 'b')
points(noise.levels, power.hhg[2,], pch = 4, type = 'b')
points(noise.levels, power.new[2,], pch = 5, type = 'b')
legend("topright",c("dcor", "MIC", "HSIC", "HHG", expression(xi[n])), pch = c(1,2,3,4,5))


plot(noise.levels, power.dcor[3,], ylim = c(0,1), xlim = c(0,1),  main = "W-shaped", xlab = "Noise Level", ylab = "Power", pch = 1, type = 'b')
points(noise.levels, power.mic[3,], pch = 2, type = 'b')
points(noise.levels, power.hsic[3,], pch = 3, type = 'b')
points(noise.levels, power.hhg[3,], pch = 4, type = 'b')
points(noise.levels, power.new[3,], pch = 5, type = 'b')
legend("topright",c("dcor", "MIC", "HSIC", "HHG", expression(xi[n])), pch = c(1,2,3,4,5))


plot(noise.levels, power.dcor[4,], ylim = c(0,1), xlim = c(0,1),  main = "Sinusoid", xlab = "Noise Level", ylab = "Power", pch = 1, type = 'b')
points(noise.levels, power.mic[4,], pch = 2, type = 'b')
points(noise.levels, power.hsic[4,], pch = 3, type = 'b')
points(noise.levels, power.hhg[4,], pch = 4, type = 'b')
points(noise.levels, power.new[4,], pch = 5, type = 'b')
legend("topright",c("dcor", "MIC", "HSIC", "HHG", expression(xi[n])), pch = c(1,2,3,4,5))

plot(noise.levels, power.dcor[5,], ylim = c(0,1), xlim = c(0,1),  main = "Circular", xlab = "Noise Level", ylab = "Power", pch = 1, type = 'b')
points(noise.levels, power.mic[5,], pch = 2, type = 'b')
points(noise.levels, power.hsic[5,], pch = 3, type = 'b')
points(noise.levels, power.hhg[5,], pch = 4, type = 'b')
points(noise.levels, power.new[5,], pch = 5, type = 'b')
legend("topright",c("dcor", "MIC", "HSIC", "HHG", expression(xi[n])), pch = c(1,2,3,4,5))


plot(noise.levels, power.dcor[6,], ylim = c(0,1), xlim = c(0,1),  main = "Heteroskedastic", xlab = "Homoskedasticity Level", ylab = "Power", pch = 1, type = 'b')
points(noise.levels, power.mic[6,], pch = 2, type = 'b')
points(noise.levels, power.hsic[6,], pch = 3, type = 'b')
points(noise.levels, power.hhg[6,], pch = 4, type = 'b')
points(noise.levels, power.new[6,], pch = 5, type = 'b')
legend("topright",c("dcor", "MIC", "HSIC", "HHG", expression(xi[n])), pch = c(1,2,3,4,5))


dev.off()