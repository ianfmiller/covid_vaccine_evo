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

## function for running analyses

do.ess.sim<-function(rUv,rLv,plot.sim=F)
{
  RE.invader.vec<-c() #build empty vector of RE values for invader strain
  for(alpha1 in A)
  {
    states<-start.states
    
    parameters0<-c(b1=b1,b2=b2,gamma=gamma,rUv=rUv,rLv=rLv,rUc=rUc,rLc=rLc,rUcv=mean(rUv,1),rLcv=mean(rLv,1),epsilon=epsilon,alpha=0.01,p=p,omega=omega,omegav=omegav)
    out0 <- ode(states, times=c(0,0), func = "derivs", parms = parameters0,dllname = "epi.model", initfunc = "initmod",nout=32,outnames=paste0("out",0:31),method = "lsoda")
    get.matricies(out0)
    Re.alpha.delta.start<-getR0(Fmat,Vmat)
    
    parameters1<-c(b1=b1,b2=b2,gamma=gamma,rUv=rUv,rLv=rLv,rUc=rUc,rLc=rLc,rUcv=mean(rUv,1),rLcv=mean(rLv,1),epsilon=epsilon,alpha=alpha1,p=p,omega=omega,omegav=omegav)
    out1 <- ode(states, times=times, func = "derivs", parms = parameters1,dllname = "epi.model", initfunc = "initmod",nout=32,outnames=paste0("out",0:31),method = "lsoda")
    if(plot.sim) {plot.simulation(out1)}
    epi.equi.states<-out1[nrow(out1),c("S","V","I_0","I_V","I_C","I_C_V","C","C_V")]
    
    out2 <- ode(epi.equi.states, times=c(0,0), func = "derivs", parms = parameters0,dllname = "epi.model", initfunc = "initmod",nout=32,outnames=paste0("out",0:31),method = "lsoda")
    get.matricies(out2)
    Re.alpha.delta.epi.equi<-getR0(Fmat,Vmat)
    
    for (alpha2 in A)
    {          
      parameters2<-c(b1=b1,b2=b2,gamma=gamma,rUv=rUv,rLv=rLv,rUc=rUc,rLc=rLc,rUcv=mean(rUv,1),rLcv=mean(rLv,1),epsilon=epsilon,alpha=alpha2,p=p,omega=omega,omegav=omegav)
      out3 <- ode(epi.equi.states, times=c(0,0), func = "derivs", parms = parameters2,dllname = "epi.model", initfunc = "initmod",nout=32,outnames=paste0("out",0:31))
      get.matricies(out3)
      RE.invader<-getR0(Fmat,Vmat)
      RE.invader.vec<-c(RE.invader.vec,RE.invader)
    }
    #print(paste0("finished alpha1 = ",alpha1))
  }
  RE.invader.mat<-matrix(RE.invader.vec,length(A),length(A),byrow = T)
  colnames(RE.invader.mat)<-A
  rownames(RE.invader.mat)<-A
  image(RE.invader.mat>=1)
  ess.result<-pip.analysis(RE.invader.mat)
  data.frame(
    "rUv"=rUv,
    "rLv"=rLv,
    "Re.alpha.delta.start"=Re.alpha.delta.start,
    "Re.alpha.delta.epi.equi"=Re.alpha.delta.epi.equi,
    "pip.motif"=paste0(na.omit(ess.result[c(1,3,5)]),collapse=" "),
    "alpha.ess"=ess.result[2],
    "repeller.ess"=ess.result[4]
  )
  print(paste0("finished rUv = ",rUv," rLv = ",rLv))
}

## set global parameters

times<-seq(0,365*400,1)
A<-seq(.0025,.2,.0005)
res<-11
rUv.steps<-rLv.steps<-seq(.5,1,length.out = res)

gamma<-1/7
epsilon<-.5
p<-50

# 10 year waning

omega<-1/(365*10)
omegav<-1/(365*10)

## alpha optim = 0.00875

### set b1, b2

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

### 90% vaccinated

start.states<-get.states(.375,.9,.001) # set startinng conditions
rUc<-.5 #convalescent class
rLc<-.75 #convalescent class

### do analysis

if(!file.exists("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega10p.vacc0.9alpha.optim0.00875.RDS"))
{
  n.cores<-detectCores()
  registerDoParallel(n.cores)
  sim.params<-data.frame("rUv"=rep(rUv.steps,each=res),"rLv"=rep(rLv.steps,times=res))
  out.data<-foreach(k = 1:nrow(sim.params), .multicombine = T, .combine = rbind, .verbose = T) %dopar% do.ess.sim(sim.params[k,"rUv"],sim.params[k,"rLv"])
  saveRDS(out.data,file="~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega10p.vacc0.9alpha.optim0.00875.RDS")
}




