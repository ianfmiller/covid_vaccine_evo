#### set up
source("covid_SVIC_functions.R")

library(deSolve)
library(rootSolve)
library(plotrix)

system("R CMD SHLIB SVIC.c")
dyn.load(paste("SVIC", .Platform$dynlib.ext, sep = ""))

#### define model parameters

gamma<-1
optim.vir.assumed<-.005*1.5 #set to either .01, .005, .0025
vir.obs<-.005*1.5
R0.assumed<-2.5*1.5
rU<-.0 #vaccinated class
rL<-0 #vaccinated class
rUn<-0 #convalescent class
rLn<-0 #convalescent class
frac_lower<-.5 # % contribution of lower respiratory infection to overall transmission

par(mfrow=c(1,3))

prop<-100
get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,100),vr=vir.obs,b2=b2,tol=1e-10)$root
x.cords<-seq(0,.015,.0001)
y.cords<-b1*x.cords^b2
plot(x.cords,y.cords,type="l",col="black",xlim=c(0,.015),xlab="virulence",ylab="transmisson rate",main="transmission rate")

prop<-66.666
get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,100),vr=vir.obs,b2=b2,tol=1e-10)$root
x.cords<-seq(0,.015,.0001)
y.cords<-b1*x.cords^b2
points(x.cords,y.cords,type="l",col="green")

prop<-10
get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,100),vr=vir.obs,b2=b2,tol=1e-10)$root
x.cords<-seq(0,.015,.0001)
y.cords<-b1*x.cords^b2
points(x.cords,y.cords,type="l",col="red")

prop<-100
x.cords<-seq(0,.015,.0001)
y.cords<-1/(x.cords*prop+gamma)
plot(x.cords,y.cords,type="l",col="black",xlab="virulence",ylab="transmisson time",main="transmission time")

prop<-66.666
x.cords<-seq(0,.015,.0001)
y.cords<-1/(x.cords*prop+gamma)
points(x.cords,y.cords,type="l",col="green")

prop<-10
x.cords<-seq(0,.015,.0001)
y.cords<-1/(x.cords*prop+gamma)
points(x.cords,y.cords,type="l",col="red")

prop<-100
get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,100),vr=vir.obs,b2=b2,tol=1e-10)$root
x.cords<-seq(0,.015,.0001)
y.cords<-sapply(x.cords,find.R0,b1=b1,b2=b2)
plot(x.cords,y.cords,type="l",col="black",xlab="virulence",ylab="fitness",main="fitness")

prop<-66.666
get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,100),vr=vir.obs,b2=b2,tol=1e-10)$root
x.cords<-seq(0,.015,.0001)
y.cords<-sapply(x.cords,find.R0,b1=b1,b2=b2)
points(x.cords,y.cords,type="l",col="green")

prop<-10
get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,100),vr=vir.obs,b2=b2,tol=1e-10)$root
x.cords<-seq(0,.015,.0001)
y.cords<-sapply(x.cords,find.R0,b1=b1,b2=b2)
points(x.cords,y.cords,type="l",col="red")

legend("bottomright",legend=c("prop = 100","prop = 66.666","prop = 10"),col=c("black","green","red"),lwd=3)
