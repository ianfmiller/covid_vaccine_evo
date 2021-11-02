# set up

## load packages,functions,  model
library(deSolve)
library(rootSolve)
library(plotrix)
library(parallel)
library(doParallel)
library(foreach)
library(doRNG)
source("~/Documents/GitHub/covid_vaccines_virulence_evolution/functions.R")
system("R CMD SHLIB ~/Documents/GitHub/covid_vaccines_virulence_evolution/epi.model.c")
dyn.load(paste("~/Documents/GitHub/covid_vaccines_virulence_evolution/epi.model", .Platform$dynlib.ext, sep = ""))

## set global parameters

times<-seq(0,365*1,1)

A<-c(seq(.0025,.025,length.out = 200))
res<-11
rUv.steps<-rLv.steps<-seq(.5,1,length.out = res)

gamma<-1/7
q=1/10
epsilon<-.5
p<-50

omega<-0
omegav<-0

optim.alpha.assumed<-.00875 #set to either .00875, .01, .02
alpha.obs<-.01
R0.assumed<-5.625

rUv<-0 #vaccinated class
rLv<-0 #vaccinated class
rUc<-0 #convalescent class
rLc<-0 #convalescent class
rUcv<-0 #vaccinated + convalescent class
rLcv<-0 #vaccinated + convalescent class

get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.alpha.assumed=optim.alpha.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,200),alpha=alpha.obs,b2=b2,tol=1e-10)$root

start.states<-get.states(.375,.9,.001) # set startinng conditions
rUc<-.5 #convalescent class
rLc<-.75 #convalescent class

# plot
par(mfrow=c(1,3),mar=c(5,5,10,5))

## pannel A
rLv<-.5
rUv<-.5

RE.invader.vec<-c() #build empty vector of RE values for invader strain
for(alpha1 in A)
{
  states<-start.states
  
  parameters0<-c(b1=b1,b2=b2,gamma=gamma,rUv=rUv,rLv=rLv,rUc=rUc,rLc=rLc,rUcv=mean(rUv,1),rLcv=mean(rLv,1),epsilon=epsilon,alpha=0.01,p=p,omega=omega,omegav=omegav,q=q)
  out0 <- ode(states, times=c(0,0), func = "derivs", parms = parameters0,dllname = "epi.model", initfunc = "initmod",nout=72,outnames=paste0("out",0:71),method = "lsoda")
  get.matricies(out0)
  Re.alpha.delta.start<-getR0(Fmat,Vmat)
  
  parameters1<-c(b1=b1,b2=b2,gamma=gamma,rUv=rUv,rLv=rLv,rUc=rUc,rLc=rLc,rUcv=mean(rUv,1),rLcv=mean(rLv,1),epsilon=epsilon,alpha=alpha1,p=p,omega=omega,omegav=omegav,q=q)
  out1 <- ode(states, times=times, func = "derivs", parms = parameters1,dllname = "epi.model", initfunc = "initmod",nout=72,outnames=paste0("out",0:71),method = "lsoda")
  epi.equi.states<-out1[nrow(out1),c("S","V","I_0","I_V","I_C","I_C_V","Q","Q_V","C","C_V")]
  
  out2 <- ode(epi.equi.states, times=c(0,0), func = "derivs", parms = parameters0,dllname = "epi.model", initfunc = "initmod",nout=72,outnames=paste0("out",0:71),method = "lsoda")
  get.matricies(out2)
  Re.alpha.delta.epi.equi<-getR0(Fmat,Vmat)
  
  for (alpha2 in A)
  {          
    parameters2<-c(b1=b1,b2=b2,gamma=gamma,rUv=rUv,rLv=rLv,rUc=rUc,rLc=rLc,rUcv=mean(rUv,1),rLcv=mean(rLv,1),epsilon=epsilon,alpha=alpha2,p=p,omega=omega,omegav=omegav,q=q)
    out3 <- ode(epi.equi.states, times=c(0,0), func = "derivs", parms = parameters2,dllname = "epi.model", initfunc = "initmod",nout=72,outnames=paste0("out",0:71))
    get.matricies(out3)
    RE.invader<-getR0(Fmat,Vmat)
    RE.invader.vec<-c(RE.invader.vec,RE.invader)
  }
  print(paste0("finished alpha1 = ",alpha1))
}
RE.invader.mat<-matrix(RE.invader.vec,length(A),length(A),byrow = T)
colnames(RE.invader.mat)<-A
rownames(RE.invader.mat)<-A
image(RE.invader.mat>=1,axes=F,col=c("grey","black"))
axis(1,at=c(0,1),labels=A[c(1,200)])
axis(2,at=c(0,1),labels=A[c(1,200)])
mtext(expression(alpha[resident]),side=1,cex=2,line=2)
mtext(expression(alpha[invader]),side=2,cex=2,line=2)
points(.46,.46,pch=1,cex=4,col="red")
mtext("A",side=3,line=1,font=2,adj=0,padj = 0)
mtext("ESS",cex=1.5,line=1)

## pannel B
rLv<-1
rUv<-.5

RE.invader.vec<-c() #build empty vector of RE values for invader strain
for(alpha1 in A)
{
  states<-start.states
  
  parameters0<-c(b1=b1,b2=b2,gamma=gamma,rUv=rUv,rLv=rLv,rUc=rUc,rLc=rLc,rUcv=mean(rUv,1),rLcv=mean(rLv,1),epsilon=epsilon,alpha=0.01,p=p,omega=omega,omegav=omegav,q=q)
  out0 <- ode(states, times=c(0,0), func = "derivs", parms = parameters0,dllname = "epi.model", initfunc = "initmod",nout=72,outnames=paste0("out",0:71),method = "lsoda")
  get.matricies(out0)
  Re.alpha.delta.start<-getR0(Fmat,Vmat)
  
  parameters1<-c(b1=b1,b2=b2,gamma=gamma,rUv=rUv,rLv=rLv,rUc=rUc,rLc=rLc,rUcv=mean(rUv,1),rLcv=mean(rLv,1),epsilon=epsilon,alpha=alpha1,p=p,omega=omega,omegav=omegav,q=q)
  out1 <- ode(states, times=times, func = "derivs", parms = parameters1,dllname = "epi.model", initfunc = "initmod",nout=72,outnames=paste0("out",0:71),method = "lsoda")
  epi.equi.states<-out1[nrow(out1),c("S","V","I_0","I_V","I_C","I_C_V","Q","Q_V","C","C_V")]
  
  out2 <- ode(epi.equi.states, times=c(0,0), func = "derivs", parms = parameters0,dllname = "epi.model", initfunc = "initmod",nout=72,outnames=paste0("out",0:71),method = "lsoda")
  get.matricies(out2)
  Re.alpha.delta.epi.equi<-getR0(Fmat,Vmat)
  
  for (alpha2 in A)
  {          
    parameters2<-c(b1=b1,b2=b2,gamma=gamma,rUv=rUv,rLv=rLv,rUc=rUc,rLc=rLc,rUcv=mean(rUv,1),rLcv=mean(rLv,1),epsilon=epsilon,alpha=alpha2,p=p,omega=omega,omegav=omegav,q=q)
    out3 <- ode(epi.equi.states, times=c(0,0), func = "derivs", parms = parameters2,dllname = "epi.model", initfunc = "initmod",nout=72,outnames=paste0("out",0:71))
    get.matricies(out3)
    RE.invader<-getR0(Fmat,Vmat)
    RE.invader.vec<-c(RE.invader.vec,RE.invader)
  }
  print(paste0("finished alpha1 = ",alpha1))
}
RE.invader.mat<-matrix(RE.invader.vec,length(A),length(A),byrow = T)
colnames(RE.invader.mat)<-A
rownames(RE.invader.mat)<-A
image(RE.invader.mat>=1,axes=F,col=c("grey","black"))
axis(1,at=c(0,1),labels=A[c(1,200)])
axis(2,at=c(0,1),labels=A[c(1,200)])
mtext(expression(alpha[resident]),side=1,cex=2,line=2)
mtext(expression(alpha[invader]),side=2,cex=2,line=2)
mtext("B",side=3,line=1,font=2,adj=0,padj = 0)
mtext("Unbounded selection for\nincreased virulence",cex=1.5,line=1)

## pannel C
rLv<-1
rUv<-1

RE.invader.vec<-c() #build empty vector of RE values for invader strain
for(alpha1 in A)
{
  states<-start.states
  
  parameters0<-c(b1=b1,b2=b2,gamma=gamma,rUv=rUv,rLv=rLv,rUc=rUc,rLc=rLc,rUcv=mean(rUv,1),rLcv=mean(rLv,1),epsilon=epsilon,alpha=0.01,p=p,omega=omega,omegav=omegav,q=q)
  out0 <- ode(states, times=c(0,0), func = "derivs", parms = parameters0,dllname = "epi.model", initfunc = "initmod",nout=72,outnames=paste0("out",0:71),method = "lsoda")
  get.matricies(out0)
  Re.alpha.delta.start<-getR0(Fmat,Vmat)
  
  parameters1<-c(b1=b1,b2=b2,gamma=gamma,rUv=rUv,rLv=rLv,rUc=rUc,rLc=rLc,rUcv=mean(rUv,1),rLcv=mean(rLv,1),epsilon=epsilon,alpha=alpha1,p=p,omega=omega,omegav=omegav,q=q)
  out1 <- ode(states, times=times, func = "derivs", parms = parameters1,dllname = "epi.model", initfunc = "initmod",nout=72,outnames=paste0("out",0:71),method = "lsoda")
  epi.equi.states<-out1[nrow(out1),c("S","V","I_0","I_V","I_C","I_C_V","Q","Q_V","C","C_V")]
  
  out2 <- ode(epi.equi.states, times=c(0,0), func = "derivs", parms = parameters0,dllname = "epi.model", initfunc = "initmod",nout=72,outnames=paste0("out",0:71),method = "lsoda")
  get.matricies(out2)
  Re.alpha.delta.epi.equi<-getR0(Fmat,Vmat)
  
  for (alpha2 in A)
  {          
    parameters2<-c(b1=b1,b2=b2,gamma=gamma,rUv=rUv,rLv=rLv,rUc=rUc,rLc=rLc,rUcv=mean(rUv,1),rLcv=mean(rLv,1),epsilon=epsilon,alpha=alpha2,p=p,omega=omega,omegav=omegav,q=q)
    out3 <- ode(epi.equi.states, times=c(0,0), func = "derivs", parms = parameters2,dllname = "epi.model", initfunc = "initmod",nout=72,outnames=paste0("out",0:71))
    get.matricies(out3)
    RE.invader<-getR0(Fmat,Vmat)
    RE.invader.vec<-c(RE.invader.vec,RE.invader)
  }
  print(paste0("finished alpha1 = ",alpha1))
}
RE.invader.mat<-matrix(RE.invader.vec,length(A),length(A),byrow = T)
colnames(RE.invader.mat)<-A
rownames(RE.invader.mat)<-A
image(RE.invader.mat>=1,axes=F,col=c("grey","black"))
axis(1,at=c(0,1),labels=A[c(1,200)])
axis(2,at=c(0,1),labels=A[c(1,200)])
mtext(expression(alpha[resident]),side=1,cex=2,line=2)
mtext(expression(alpha[invader]),side=2,cex=2,line=2)
mtext("C",side=3,line=1,font=2,adj=0,padj = 0)
mtext("Global eradication",cex=1.5,line=1)
legend("topright",legend=c("resident strain wins","invader strain wins"),pch=15,col=c("grey","black"),cex=2)

#save to image with dimensions 1289x492

