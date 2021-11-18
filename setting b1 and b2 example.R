# set up
source("functions.R")

library(deSolve)
library(rootSolve)
library(plotrix)

system("R CMD SHLIB epi.model.c")
dyn.load(paste("epi.model", .Platform$dynlib.ext, sep = ""))

# define model parameters fo rsetting b1 and b2

rUv<-.0 #vaccinated class
rLv<-0 #vaccinated class
rUc<-0 #convalescent class
rLc<-0 #convalescent class
rUcv<-0 #vaccinated + convalescent class
rLcv<-0 #vaccinated + convalescent class

gamma<-1/7
epsilon<-.5
p<-50
omega<-0
omegav<-0
mu<-1/(73*365)
f<-0

optim.alpha.assumed<-.01 #set to either .00875, .01, .02
alpha.obs<-.01

R0.assumed<-5.625


#set trade-off (b1 + b2) according to assumed ES virulence, observed virulence, and R0 at observed virulence

## set trade-off shape (b2 only) according to assumed ES virulence

get.states(0,0,0,0)

b2<-uniroot(b2.search,c(0,1),b1=1,optim.alpha.assumed=optim.alpha.assumed,tol=1e-15)$root

## set trade-off scaling (b1) so that R0=2.5 at vir.obs

b1<-uniroot(R0.search,c(0,200),vr=alpha.obs,b2=b2,tol=1e-10)$root

# checks for setting b1 and b2

##check to ensure that fit b2 is effectively independent of b1

uniroot(b2.search,c(0,1),b1=1,optim.alpha.assumed=optim.alpha.assumed,tol=1e-10)$root
uniroot(b2.search,c(0,1),b1=2,optim.alpha.assumed=optim.alpha.assumed,tol=1e-10)$root
uniroot(b2.search,c(0,1),b1=3,optim.alpha.assumed=optim.alpha.assumed,tol=1e-10)$root

#check to ensure that b1 and b2 give desired R0 at vir.obs
find.R0(alpha.obs,b1,b2)

#check that optim.vir.assumed is in fact optimal given b1 and b2, and that b1 and b2 and give desired R0 at vir.obs

vs<-seq(.005,.05,.00001)
rs<-unlist(lapply(vs,find.R0,b1=b1,b2=b2))
plot(vs,rs,ylim=c(5.62,5.63),xlab=expression(alpha),ylab=expression(R[0]))
abline(v=optim.alpha.assumed)
abline(h=max(rs))
abline(v=alpha.obs)
abline(h=R0.assumed)