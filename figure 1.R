#### set up

source("covid_SVIC_functions.R")

library(deSolve)
library(rootSolve)
library(viridis)

system("R CMD SHLIB SVIC.c")
dyn.load(paste("SVIC", .Platform$dynlib.ext, sep = ""))

#### define model parameters

gamma<-1

rU<-.0 #vaccinated class
rL<-0 #vaccinated class
rUn<-0 #convalescent class
rLn<-0 #convalescent class

frac_lower<-.5 # % contribution of lower respiratory infection to overall transmission

R0.assumed<-2.5

vir.obs<-.005

### analysis
optim.vir.assumed<-.0025

get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,100),vr=vir.obs,b2=b2,tol=1e-10)$root

get.states(0,0,0)
frac_lower<-.5
vs<-seq(0,.015,.00001)

rU<-0
rL<-0
trans.rates1<-(frac_lower)*b1*(vs*(1-rL))^b2+(1-frac_lower)*b1*(vs*(1-rU))^b2
death.rates1<-vs*(1-rL)
out<-find.optim.vir(vs,b1,b2,rU,rL,0,0)
#optim.vir1<-out[[1]]
R0s1<-out[[2]]
#index1<-which(vs==optim.vir1)


optim.vir.assumed<-.005 

get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,100),vr=vir.obs,b2=b2,tol=1e-10)$root

get.states(0,0,0)
frac_lower<-.5
vs<-seq(0,.015,.00001)

rU<-0
rL<-0
trans.rates2<-(frac_lower)*b1*(vs*(1-rL))^b2+(1-frac_lower)*b1*(vs*(1-rU))^b2
death.rates2<-vs*(1-rL)
out<-find.optim.vir(vs,b1,b2,rU,rL,0,0)
#optim.vir2<-out[[1]]
R0s2<-out[[2]]
#index2<-which(vs==optim.vir2)

optim.vir.assumed<-.01 

get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,100),vr=vir.obs,b2=b2,tol=1e-10)$root

get.states(0,0,0)
frac_lower<-.5
vs<-seq(0,.015,.00001)

rU<-0
rL<-0
trans.rates3<-(frac_lower)*b1*(vs*(1-rL))^b2+(1-frac_lower)*b1*(vs*(1-rU))^b2
death.rates3<-vs*(1-rL)
out<-find.optim.vir(vs,b1,b2,rU,rL,0,0)
#optim.vir3<-out[[1]]
R0s3<-out[[2]]
#index3<-which(vs==optim.vir3)

par(mfrow=c(1,3),mar=c(5,5,4,1))
plot(vs,trans.rates1,col=viridis(3)[1],ylim=c(2.49,2.54),xlab="virulence",ylab="transmisison rate",type="l",lwd=6,cex.lab=2,cex.axis=1.2)
points(vs,trans.rates2,col=viridis(3)[2],type="l",lwd=6)
points(vs,trans.rates3,col=viridis(3)[3],type="l",lwd=6)
mtext("A",side=3,line=1,font=2,adj=0,padj = 0)

plot(vs,1/(gamma+death.rates1),col=viridis(3)[1],ylim=c(.98,1),xlab="virulence",ylab="transmisison time",type="l",lwd=6,lty=1,lend=1,cex.lab=2,cex.axis=1.2)
points(vs,1/(gamma+death.rates2),col=viridis(3)[2],type="l",lwd=6,lty=c("21"),lend=1)
points(vs,1/(gamma+death.rates3),col=viridis(3)[3],type="l",lwd=6,lty=c("12"),lend=1)
mtext("B",side=3,line=1,font=2,adj=0,padj = 0)

plot(vs,R0s1,col=viridis(3)[1],ylim=c(2.495,2.505),xlab="virulence",ylab="fitness",type="l",lwd=6,cex.lab=2,cex.axis=1.2)
points(vs,R0s2,col=viridis(3)[2],type="l",lwd=6)
points(vs,R0s3,col=viridis(3)[3],type="l",lwd=6)
mtext("C",side=3,line=1,font=2,adj=0,padj = 0)
