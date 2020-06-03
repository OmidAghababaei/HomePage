########################################################################
############# Prior and Posterior Beta distribution, page 6 
########################################################################
X11(w=8,h=5); par(mfrow=c(1,3));
tv=seq(0,1,0.01)
plot(c(0,1),c(0,1.2),type="n",xlab="theta",ylab="prior density",main="Beta(1,1)")
lines(c(0,1),c(1,1),lty=1,lwd=3)
plot(c(0,1),c(0,1.7),type="n",xlab="theta",ylab="prior density",main="Beta(2,2)")
lines(tv,3*2*tv*(1-tv),lty=1,lwd=3)
plot(c(0,1),c(0,2.7),type="n",xlab="theta",ylab="prior density",main="Beta(5,5)")
lines(tv,630*(tv^4)*(1-tv)^4,lty=1,lwd=3)
########################################################################
############# Prior and Posterior Bernoulli Model, page 11
########################################################################
X11(w=8,h=5); par(mfrow=c(1,1)); 
plot(c(0,1),c(0,3),type="n",xlab="theta",ylab="density")
lines(c(0,1),c(1,1),lty=1,lwd=3); tv=seq(0,1,0.01) 
lines(tv,3*(1-tv)^2,lty=2,lwd=3,col="red") 
lines(tv,3*2*tv*(1-tv),lty=3,lwd=3,col="blue") 
lines(tv,3*tv^2,lty=4,lwd=3,col="green")
legend(0.3,3,legend=c("prior","posterior if sum(y_i)=0","posterior if sum(y_i)=1", 
    "posterior if sum(y_i)=2"), col=c("black","red","blue","green"),lty=c(1,2,3,4),lwd=rep(2,4))
########################################################################
############# Prior and Posterior normal Model, page 15 
########################################################################
X11(w=8,h=5); par(mfrow=c(1,1)); mu0=5; tau0=0.5; sig0=1
y = c(8.4, 10.1, 9.4); n = length(y);  
k=n*tau0^2/(n*tau0^2+sig0^2); k # 0.4285714 
ybar=mean(y); ybar # 9.3
mus = (1-k)*mu0 + k*ybar;  sigs2=k*sig0^2/n
c(mus,sigs2) # 6.8428571 0.1428571
muv=seq(0,15,0.01);prior = dnorm(muv,mu0,tau0); 
post=dnorm(muv,mus,sqrt(sigs2)) 
like = dnorm(muv,ybar,sig0/sqrt(n)) 
HPD=mus+c(-1,1)*qnorm(0.975)*sqrt(sigs2)
HPD # 6.102060 7.583654
plot(c(0,11),c(-0.1,1.3),type="n",xlab="",ylab="density/likelihood") 
lines(muv,prior,lty=1,lwd=2); lines(muv,like,lty=2,lwd=2) 
lines(muv,post,lty=3,lwd=2) 
points(c(mu0,ybar,mus),c(0,0,0),pch=c(1,2,4),cex=rep(1.5,3),lwd=2) 
points(HPD,c(0,0),pch=rep(16,2),cex=rep(1.5,2))
legend(0,1.3,c("Prior density","Likelihood function (normalized)",
               "Posterior density"),lty=c(1,2,3),lwd=c(2,2,2))
legend(0,0.7,c("Prior mean","Sample mean (MLE)","Posterior mean",
               "95% HPD Interval"), pch=c(1,2,4,16),pt.cex=rep(1.5,4),pt.lwd=rep(2,4)) 
text(10.8,-0.075,"m", vfont=c("serif symbol","italic"), cex=1.5)
########################################################################
############# Example 1: the N-R algorithm , page 39 
########################################################################
NR <- function(th,J=5){
  thvec <- th; for(j in 1:J){
    num <- th^6-1/2;den <- 6*th^5
    th <- th - num/den; thvec <- c(thvec,th) }; thvec } 
# theta’s posterior cdf minus 1/2 (numerator) 
# theta’s posterior pdf (denominator)
options(digits=4)
NR(th=1,J=6)  
NR(th=0.8,J=6)  
0.8909-(0.8909^6-0.5)/(6*0.8909^5) # 0.8909 
########################################################################
############# Example 2: the N-R algorithm , page 40 
########################################################################
options(digits=6); t=0; tv=t; for(j in 1:7){ t=t-(t^2-exp(t))/(2*t-exp(t)) 
tv=c(tv,t) }; tv; t^2-exp(t)  
(-0.703467)^2-exp(-0.703467) # -8.03508e-07
t=1; tv=t; for(j in 1:7){ t=t-(t^2-exp(t))/(2*t-exp(t)); tv=c(tv,t) }; tv
t=-5; tv=t; for(j in 1:7){ t=t-(t^2-exp(t))/(2*t-exp(t)); tv=c(tv,t) }; tv
tvec=seq(-6,2,0.01); gvec= tvec^2-exp(tvec) 
X11(w=8,h=4.5); par(mfrow=c(1,1)) 
plot(tvec,gvec,type="l",lwd=2,xlab="t",ylab="g(t)", main="") 
abline(h=0,v=t); points(tv, tv^2-exp(tv),pch=16)
text( tv[1:4], tv[1:4]^2-exp(tv[1:4]) + 3, 0:3)
########################################################################
############# Example 3: the N-R algorithm , page 43 
########################################################################
cat(" Intial values= {2/3, 0.5, 0.9, 0.1, 0.614272}; Repeat= {7,20,1}")
options(digits=6); p=2/3; pv=p; for(j in 1:7){
p = p - (4*p^3-3*p^4-1/2)/(12*p^2-12*p^3); pv=c(pv,p) };round(pv,3) 

p=0.5; pv=p; for(j in 1:7){ p = p - (4*p^3-3*p^4-1/2)/(12*p^2-12*p^3); pv=c(pv,p) };  
round(pv,3) 

p=0.9; pv=p; for(j in 1:7){ p = p - (4*p^3-3*p^4-1/2)/(12*p^2-12*p^3); pv=c(pv,p) }; 
round(pv,3) 

p=0.1; pv=p; for(j in 1:7){ p = p - (4*p^3-3*p^4-1/2)/(12*p^2-12*p^3); pv=c(pv,p) };  
round(pv,3)  

p=0.1; pv=p; for(j in 1:20){ p = p - (4*p^3-3*p^4-1/2)/(12*p^2-12*p^3); pv=c(pv,p) };  
round(pv,3) 

4*(0.614272)^3-3*(0.614272)^4 # 0.499999 
########################################################################
############# Example 3: the N-R algorithm Plot, page 45 
########################################################################
pvec=seq(-0.5,1.4,0.005); 
Fvec = 4*pvec^3-3*pvec^4 
Fvec[pvec<=0] = 0; 
Fvec[pvec>=1] = 1

X11(w=8,h=4.5); par(mfrow=c(2,1))

plot(pvec,Fvec,type="l",lwd=3,xlab="p",ylab="F(p|x)", main="Posterior cdf and median of p") 
abline(h=0.5,v=0.614272,lty=3); points(0.614272,0.5,pch=16, cex=1.2) 
abline(h=c(0,1),lty=3); abline(v=c(0,1),lty=3) 
gvecwrong=4*pvec^3-3*pvec^4-0.5

plot(pvec, gvecwrong,type="n",lwd=2,xlab="p",ylab="g(p) = F(p|x) - 1/2", 
     main="Posterior median of p and the other root of g")
lines(pvec,Fvec-0.5,lwd=3)
lines(pvec[pvec<0], gvecwrong[pvec<0],lty=2,lwd=3) 
lines(pvec[pvec>1], gvecwrong[pvec>1],lty=2,lwd=3) 
abline(v=c(0.614272, 1.24748),lty=3); abline(h=0,lty=3) 
points(c(0.614272, 1.24748),c(0,0),pch=16,cex=1.2) 
abline(h=c(-0.5,0,0.5),lty=3); abline(v=c(0,1),lty=3)
################################################################################
############# Example : Finding a HPD via the multivariate N-R algorithm, page 48 
#################################################################################
cat(" Intial Values: 0.5,0.3 and Repeat=7")
gfun = function(a,b){
  g1=pgamma(b,2,1)-pgamma(a,2,1)-0.8; g2=dgamma(b,2,1)-dgamma(a,2,1); c(g1,g2) }
gpfun = function(a,b){ m11=-dgamma(a,2,1); m12=dgamma(b,2,1) 
m21=exp(-a)*(a-1); m22=exp(-b)*(1-b) 
matrix(c(m11,m12,m21,m22),nrow=2,byrow=T) }
gvec=c(0.5,3); gmat=gvec; for(j in 1:7){ 
a=gvec[1]; b=gvec[2]
gvec = gvec - solve(gpfun(a,b)) %*% gfun(a,b) 
gmat = cbind(gmat,gvec) }
options(digits=4); gmat
# Checks:
cat(" Values of a and b and the Corresponding value on the Gamma curves (c)")
c(a,b,dgamma(c(a,b),2,1)) # 0.167300 3.080291 0.141527 0.141527 
cat(" Areas under the gamma curve before and between  a and  b")
c(pgamma(a,2,1), pgamma(b,2,1), pgamma(b,2,1) - pgamma(a,2,1)) 
################################################################################
############# Example : Finding a HPD via the multivariate N-R algorithm Plot, page 49 
#################################################################################
gfun = function(a,b){
  g1=pgamma(b,2,1)-pgamma(a,2,1)-0.8; g2=dgamma(b,2,1)-dgamma(a,2,1); c(g1,g2) }
gpfun = function(a,b){ m11=-dgamma(a,2,1); m12=dgamma(b,2,1) 
m21=exp(-a)*(a-1); m22=exp(-b)*(1-b) 
matrix(c(m11,m12,m21,m22),nrow=2,byrow=T) }
gvec=c(0.5,3); gmat=gvec; for(j in 1:7){ 
a=gvec[1]; b=gvec[2]
gvec = gvec - solve(gpfun(a,b)) %*% gfun(a,b) 
gmat = cbind(gmat,gvec) }
options(digits=6);  
lamv=seq(0,5,0.01); fv=dgamma(lamv,2,1)
X11(w=8,h=3.5); par(mfrow=c(1,1)) 
plot(lamv,fv,type="l",lwd=3,xlab="lambda",ylab="f(lambda|x)", main=" ") 
abline(h=c(dgamma(a,2,1)),v=c(a,b),lty=1)
################################################################################
############# Example : Integration techniques, page 52 
#################################################################################
INTEG <- function(xvec, yvec, a = min(xvec), b = max(xvec)){
  fit <- smooth.spline(xvec, yvec)
  spline.f <- function(x){predict(fit, x)$y } 
  integrate(spline.f, a, b)$value }
mu=8; sig=3; c = 10; options(digits=6)
PXpos = (1-pnorm((c-mu)/sig)) 
gfun=function(x){ x * dnorm(x,mu,sig) / PXpos } 
integrate(gfun,c,20)$value   
integrate(gfun,c,30)$value  
xvec <- seq(c,20,0.1); gvec <- gfun(xvec); INTEG(xvec,gvec,c,20)  
xvec <- seq(c,30,0.1); gvec <- gfun(xvec); INTEG(xvec,gvec,c,30)  
################################################################################
#############Example : Specification of Parameters using optim() function, page 55 
#################################################################################
options(warn=-1)
options(digits=5); a=0.5; b=1; alp=0.05; 
fun=function(v,alp=0.05,a=0.5,b=1){(pgamma(1/a^2,v[1],v[2])-(1-alp/2))^2+(pgamma(1/b^2,v[1],v[2])-(alp/2))^2 }
res0=optim(par=c(0.2,6),fn=fun)$par
res0 # 8.4764 3.7679
pgamma(c(1/b^2,1/a^2),res0[1],res0[2]) # 0.025048 0.975104 Close
res=optim(par=res0,fn=fun)$par; res # 8.4748 3.7654 
pgamma(c(1/b^2,1/a^2),res[1],res[2]) # 0.025 0.975 Correct
res2=optim(par=c(6,3),fn=fun)$par; res2 # 8.4753 3.7655 
pgamma(c(1/b^2,1/a^2),res2[1],res2[2]) # 0.024992 0.974996 Close
res3=optim(par=res2,fn=fun)$par; res3 # 8.4748 3.7654 
pgamma(c(1/b^2,1/a^2),res3[1],res3[2]) # 0.025 0.975 Correct  
# Check areas under the last curve
func=function(t){ dgamma(1/t^2,res[1],res[2])*2/t^3 } 
integrate(func,lower=0,upper=Inf)$value # 1 Correct 
integrate(func,lower=0,upper=0.5)$value # 0.025 Correct 
integrate(func,lower=1,upper=Inf)$value # 0.025 Correct
################################################################################
#############Example : Specification of Parameters using optim() function Plot, page 57 
#################################################################################
options(digits=5); a=0.5; b=1; alp=0.05; 
par(mfrow=c(3,1)); tv=seq(0,10,0.01)
plot(tv, dgamma(tv,res[1],res[2]),type="l",lwd=2, xlim=c(0,6), xlab="lambda",ylab="density"); 
abline(v=c(1/a^2,1/b^2));
abline(h=0,lty=3)
plot(tv,dgamma(1/tv,res[1],res[2])/tv^2, type="l", lwd=2, xlim=c(0,1.5), xlab="sigma^2",ylab="density");
abline(v=c(a^2,b^2)); abline(h=0,lty=3)
plot(tv,dgamma(1/tv^2,res[1],res[2])*2/tv^3, type="l", lwd=2, xlim=c(0.35,1.4), xlab="sigma",ylab="density"); 
abline(v=c(a,b)); abline(h=0,lty=3)
################################################################################
#############Example : M.C inference - Normal model, page 63
#################################################################################
mu0=2.75; tau0=0.5; sig0=1;y=c(2.1, 3.2, 5.2, 1.7);n=length(y); 
ybar=mean(y);k=n*tau0^2/(n*tau0^2+sig0^2/n); k ;ybar=mean(y); ybar 
mus = (1-k)*mu0 + k*ybar; sigs2=k*sig0^2/n;c(mus,sigs2)
J=1000; set.seed(144); zj<-rnorm(J);muj=ybar+sqrt(sigs2)*zj
mubar=mean(muj); muCI=mubar + c(-1,1)*qnorm(0.975)*sd(muj)/sqrt(J) 
muHPD=quantile(muj,c(0.025,0.975));c(mubar,muCI,muHPD) 
muhat=mus; muHPDtrue= mus+ sqrt(sigs2)*qnorm(c(0.025,0.975)) 
c(muhat,muHPDtrue) 
################################################################################
#############Example : M.C inference - Normal model Plot, page 64
#################################################################################
mu0=2.75; tau0=0.5; sig0=1;y=c(2.1, 3.2, 5.2, 1.7);n=length(y);s=sd(y)
ybar=mean(y);k=n*tau0^2/(n*tau0^2+sig0^2/n);ybar=mean(y);
mus = (1-k)*mu0 + k*ybar; sigs2=k*sig0^2/n;#c(mus,sigs2)
J=1000; set.seed(144); zj<-rnorm(J);muj=ybar+sqrt(sigs2)*zj
mubar=mean(muj); muCI=mubar + c(-1,1)*qnorm(0.975)*sd(muj)/sqrt(J) 
muHPD=quantile(muj,c(0.025,0.975));#c(mubar,muCI,muHPD) 
muhat=mus; muHPDtrue= mus+ sqrt(sigs2)*qnorm(c(0.025,0.975)) 
#c(muhat,muHPDtrue) 
X11(w=8,h=5); par(mfrow=c(1,1))
hist(muj,prob=T,xlab="mu",xlim=c(1,4.5), ylim=c(0,1),main="", breaks=seq(-20,20,0.25))
muvec=seq(-20,-20,0.01);
postvec=dnorm((muvec-mus)/sqrt(sigs2))/(s/sqrt(n))
lines(muvec,postvec, lty=2,lwd=3)
lines(density(muj),lty=1,lwd=3) 
abline(v=c(mubar,muCI,muHPD),lty=2,lwd=3)
abline(v=c(muhat, muHPDtrue) , lty=1,lwd=3)
legend(1,1,c("Monte Carlo estimates","Exact posterior estimates"),lty=c(2,1),lwd=c(3,3),bg="white")
################################################################################
############# Example : Monte Carlo prediction - binomial-beta model, page 68
#################################################################################
options(digits=5)
n=50; y=32; alp=1;bet=1; a=alp+y; b=bet+n-y; m=10; J=10000 
set.seed(443); tv=rbeta(J,a,b); xv=rbinom(J,m,tv) 
phat=length(xv[xv>=6])/J; CI=phat+c(-1,1)*qnorm(0.975)*sqrt(phat*(1-phat)/J)
c(phat,CI) 
xvec=0:m; 
fxgiveny= choose(m,xvec)*beta(y+xvec+alp,n-y+m-xvec+bet)/beta(y+alp,n-y+bet)
sum(fxgiveny)  # Check it is a proper probability function
sum(fxgiveny[xvec>=6]) # Check the exact value of P(X >= 6)
################################################################################
############# Example : Empirical Bayes - Bernoulli Model, page 77
#################################################################################
options(warn=-1)
options(digits=5); n=20;ybar=5/20;lamb=seq(0,20,0.01)
m_lamb<-gamma(2*lamb)*gamma(n*ybar+lamb)*gamma(n*(1-ybar)+lamb)/(((gamma(lamb))^2)*gamma(n+2*lamb))
plot(lamb,m_lamb,type="l",lwd=2, xlim=c(0,20), xlab="lambda",ylab="prior predictive")
abline(v=2.3,lty=2,col="red",lwd=4)
###############################