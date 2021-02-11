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

R0.assumed<-3.75

vir.obs<-.0075

### analysis
#vs<-seq(0.0055,.017,.00001)
vs<-seq(.0025,.017,.00001)
vs.alt<-seq(0,40,.01)
#### scenario 1
optim.vir.assumed<-.0075
prop<-66.666

get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,200),vr=vir.obs,b2=b2,tol=1e-10)$root

trans.rates1<-b1*(vs-.0025)^b2
trans.rates1.alt<-b1*(vs.alt-.0025)^b2
death.rates1<-1/(vs*prop+gamma)
out<-find.optim.vir(vs,b1,b2,rU,rL,0,0)
R0s1<-out[[2]]




#### scenario 2
optim.vir.assumed<-.0075*1.5
prop<-66.666

get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,200),vr=vir.obs,b2=b2,tol=1e-10)$root

trans.rates2<-b1*(vs-.0025)^b2
trans.rates2.alt<-b1*(vs.alt-.0025)^b2
death.rates2<-1/(vs*prop+gamma)
out<-find.optim.vir(vs,b1,b2,rU,rL,0,0)
R0s2<-out[[2]]

#### scenario 3
optim.vir.assumed<-.0075*2
prop<-66.666

get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,200),vr=vir.obs,b2=b2,tol=1e-10)$root

trans.rates3<-b1*(vs-.0025)^b2
trans.rates3.alt<-b1*(vs.alt-.0025)^b2
death.rates3<-1/(vs*prop+gamma)
out<-find.optim.vir(vs,b1,b2,rU,rL,0,0)
R0s3<-out[[2]]



### plot

par(mfrow=c(1,3),mar=c(5,5,4,1))
plot(vs,trans.rates1,xlim=c(0.0055,.017),ylim=c(1.75,5),xlab="virulence",ylab="transmisison rate",type="n",cex.lab=2,cex.axis=1.2)
abline(v=.0075,col=viridis(3,alpha=.5)[1],lwd=4,lty=2)
abline(v=.01125,col=viridis(3,alpha=.5)[2],lwd=4,lty=2)
abline(v=.015,col=viridis(3,alpha=.5)[3],lwd=4,lty=2)
points(vs,trans.rates1,col=viridis(3)[1],type="l",lwd=12,lty=1)
points(vs,trans.rates2,col=viridis(3)[2],type="l",lwd=12,lty=1)
points(vs,trans.rates3,col=viridis(3)[3],type="l",lwd=12,lty=1)
mtext("A",side=3,line=1,font=2,adj=0,padj = 0)

#par(fig = c(1/9,1/3-.001, 0.05, .43), new = T) 
#plot(vs.alt,trans.rates1.alt,ylim=c(0,4000),xlab="",ylab="",type="n",cex.lab=2,cex.axis=1.2)
#rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
#points(vs.alt,trans.rates1.alt,col=viridis(3)[1],type="l",lwd=12,lty=1)
#points(vs.alt,trans.rates2.alt,col=viridis(3)[2],type="l",lwd=12,lty=1)
#points(vs.alt,trans.rates3.alt,col=viridis(3)[3],type="l",lwd=12,lty=1)

#par(fig = c(1/3,2/3, 0, 1), new = T)
plot(vs,death.rates1,xlim=c(0.0055,.017),ylim=c(.5,2),xlab="virulence",ylab="transmisison time",type="n",cex.lab=2,cex.axis=1.2)
abline(v=.0075,col=viridis(3,alpha=.5)[1],lwd=4,lty=2)
abline(v=.01125,col=viridis(3,alpha=.5)[2],lwd=4,lty=2)
abline(v=.015,col=viridis(3,alpha=.5)[3],lwd=4,lty=2)
points(vs,death.rates1,col=viridis(3)[1],type="l",lwd=12,lty=1)
points(vs,death.rates2,col=viridis(3)[2],type="l",lwd=6,lty=1)
points(vs,death.rates3,col=viridis(3)[3],type="l",lwd=4,lty=2)
mtext("B",side=3,line=1,font=2,adj=0,padj = 0)

#par(fig = c(2/3,3/3, 0, 1), new = T)
plot(vs,R0s1,xlim=c(0.0055,.017),ylim=c(3.4,4.5),xlab="virulence",ylab="fitness",type="n",cex.lab=2,cex.axis=1.2)
abline(v=.0075,col=viridis(3,alpha=.5)[1],lwd=4,lty=2)
abline(v=.01125,col=viridis(3,alpha=.5)[2],lwd=4,lty=2)
abline(v=.015,col=viridis(3,alpha=.5)[3],lwd=4,lty=2)
points(vs,R0s1,col=viridis(3)[1],type="l",lwd=12)
points(vs,R0s2,col=viridis(3)[2],type="l",lwd=12)
points(vs,R0s3,col=viridis(3)[3],type="l",lwd=12)
mtext("C",side=3,line=1,font=2,adj=0,padj = 0)

legend("bottomright",legend=c("optim vir = 0.0075","optim vir = 0.01125","optim vir = 0.015"),col=viridis(3),lwd=12,cex=1.5)


