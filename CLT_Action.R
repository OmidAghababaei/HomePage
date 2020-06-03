remove(list=ls())
par(mfrow=c(3,2))
#######################
nsim <- 10^3
#n <- 100 #12
for(try in 1:6){
  
  if(try==1){n<-10}
  if(try==2){n<-50}
  if(try==3){n<-100}
  if(try==4){n<-200}
  if(try==5){n<-400}
  if(try==6){n<-1000}
  

x <- matrix(runif(n*nsim), nrow=nsim, ncol=n)
xbar <- rowMeans(x)
#hist(xbar)

walk <- rnorm(n, 4, 1) ; bus <- runif(n, 4, 16)
ride <- rpois(n, 8); climb <- rgamma(n, shape = 6, scale = 0.5)
fall <- rexp(n, rate = 4)

DT <- data.frame(walk, bus, ride, climb, fall)
DT$transit_time <- apply(DT,1,sum)

mean.walk=4; mean.bus=(4+16)/2; mean.ride=8; mean.climb=6*0.5; mean.fall=1/4
var.walk=1; var.bus=(16-4)^2/12; var.ride=8; var.climb=6*0.5^2; var.fall=1/4^2

mean.transit <- mean.walk+mean.bus+mean.ride+mean.climb+mean.fall
var.transit <- var.walk+var.bus+var.ride+var.climb+var.fall


hist(DT$transit_time, freq=FALSE,main = "")
curve(dnorm(x,mean=mean.transit,sd=sqrt(var.transit)),0,50,add=TRUE,lwd=2,col="red")
title(paste("Sample  Size =",n))


}



