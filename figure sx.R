#### set up

source("covid_SVIC_functions.R")

library(deSolve)
library(rootSolve)
library(viridis)

system("R CMD SHLIB SVIC.c")
dyn.load(paste("SVIC", .Platform$dynlib.ext, sep = ""))

#### define model parameters

gamma<-1/7

rU<-.0 #vaccinated class
rL<-0 #vaccinated class
rUn<-0 #convalescent class
rLn<-0 #convalescent class

frac_lower<-.5 # % contribution of lower respiratory infection to overall transmission

R0.assumed<-2.5*1.5

vir.obs<-.005*1.5

### analysis
vs<-seq(0,.02,.00001)
#### scenario 1
optim.vir.assumed<-.005*1.5
prop<-66.666

get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,200),vr=vir.obs,b2=b2,tol=1e-10)$root

prop<-10

trans.rates1<-b1*vs^b2
death.rates1<-1/(vs*prop+gamma)
out<-find.optim.vir(vs,b1,b2,rU,rL,0,0)
R0s1<-out[[2]]

#### scenario 2
optim.vir.assumed<-.005*1.5
prop<-66.666

get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,200),vr=vir.obs,b2=b2,tol=1e-10)$root

prop<-66.666

trans.rates2<-b1*vs^b2
death.rates2<-1/(vs*prop+gamma)
out<-find.optim.vir(vs,b1,b2,rU,rL,0,0)
R0s2<-out[[2]]

#### scenario 3
optim.vir.assumed<-.005*1.5
prop<-66.666

get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,200),vr=vir.obs,b2=b2,tol=1e-10)$root

prop<-100

trans.rates3<-b1*vs^b2
death.rates3<-1/(vs*prop+gamma)
out<-find.optim.vir(vs,b1,b2,rU,rL,0,0)
R0s3<-out[[2]]

### plot

par(mfrow=c(1,3),mar=c(5,5,4,1))
plot(vs,trans.rates1,col=viridis(3)[1],ylim=c(0,6),xlab="virulence",ylab="transmisison rate",type="l",lwd=12,cex.lab=2,cex.axis=1.2)
points(vs,trans.rates2,col=viridis(3)[2],type="l",lwd=4,lty=5)
points(vs,trans.rates3,col=viridis(3)[3],type="l",lwd=4,lty=3)
mtext("A",side=3,line=1,font=2,adj=0,padj = 0)
legend("topleft",legend=c("prop = 10","prop = 66.666","prop = 100"),col=viridis(3),lwd=12,cex=2)


plot(vs,death.rates1,col=viridis(3)[1],ylim=c(0,7),xlab="virulence",ylab="transmisison time",type="l",lwd=12,lty=1,cex.lab=2,cex.axis=1.2)
points(vs,death.rates2,col=viridis(3)[2],type="l",lwd=12,lty=1)
points(vs,death.rates3,col=viridis(3)[3],type="l",lwd=12,lty=1)
mtext("B",side=3,line=1,font=2,adj=0,padj = 0)

plot(vs,R0s1,col=viridis(3)[1],ylim=c(0,15),xlab="virulence",ylab="fitness",type="l",lwd=12,cex.lab=2,cex.axis=1.2)
points(vs,R0s2,col=viridis(3)[2],type="l",lwd=12)
points(vs,R0s3,col=viridis(3)[3],type="l",lwd=12)
mtext("C",side=3,line=1,font=2,adj=0,padj = 0)

