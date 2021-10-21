# set up

## load packages,functions,  model
library(deSolve)
library(rootSolve)
library(plotrix)
source("covid_SVIC_functions.R")
system("R CMD SHLIB SVIC.c")
dyn.load(paste("SVIC", .Platform$dynlib.ext, sep = ""))

## plotting tools

col.trans<-.6

s.blues1<-hsv(.666,1,1,seq(1,.001,length.out = 200)^col.trans)
s.reds1<-hsv(1,1,1,seq(.001,1,length.out = 200)^col.trans)
s.colors1<-c(s.blues1,"white",s.reds1)
s.col.vals1<-seq(-7.5,5.5,length.out = 401)

s.blues2<-hsv(.666,1,1,seq(1,.001,length.out = 200)^col.trans)
s.reds2<-hsv(1,1,1,seq(.001,1,length.out = 200)^col.trans)
s.colors2<-c(s.blues2,"white",s.reds2)
s.col.vals2<-seq(-7.5,7.5,length.out = 401)

s.blues3<-hsv(.666,1,1,seq(1,.001,length.out = 200)^col.trans)
s.reds3<-hsv(1,1,1,seq(.001,1,length.out = 200)^col.trans)
s.colors3<-c(s.blues3,"white",s.reds3)
s.col.vals3<-seq(-7.5,7.5,length.out = 401)

s.colors<-list(s.colors1,s.colors2,s.colors3)
s.col.vals<-list(s.col.vals1,s.col.vals2,s.col.vals3)

res<-21

## set plot window

#layout(matrix(c(1,1,2,2,3,3,10,4,4,5,5,6,6,10,7,7,8,8,9,9,10),3,7,byrow = T))
#par(mar=c(2,2,2,2),oma=c(4,8,4,0))

# analysis

## set global parameters

color.index<-1 #index to match colors and color values to right analysis

virulence.steps<-seq(.003,.35,.005)
times<-seq(0,365*100,1)
         
gamma<-1/7
epsilon<-.5
p<-50
#mu<-1/(73*365)
mu<-0

## 1st scenario

### set scenario parameters

omega<-1/(10*365)
omegav<-1/(10*365)
f<-1

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

get.states(0,0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.alpha.assumed=optim.alpha.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,200),alpha=alpha.obs,b2=b2,tol=1e-10)$root

### set scenario parameters

rUc<-.8 #convalescent class
rLc<-.95 #convalescent class
rUcv<-.99 #vaccinated + convalescent class
rLcv<-.99 #vaccinated + convalescent class
start.states<-get.states(.25,.1,.01,.5) # set startinng conditions

if(!file.exists("~/Documents/GitHub/covid_vaccine_evo/sim.data/rUc0.5rLc0.75p.vacc0.1alpha.optim0.00875.RDS"))
{
  virulence.steps.mod<-seq(.05,.09,.001)
  data<-list()
  
  index<-1
  
  for (rUv in seq(0,1,length.out = res))
  {
    for (rLv in seq(0,1,length.out = res))
    {
      RE.invader.vec<-c() #build empty vector of RE values for invader strain
      for(alpha1 in virulence.steps.mod)
      {
        states<-start.states
        
        parameters0<-c(b1=b1,b2=b2,gamma=gamma,rUv=rUv,rLv=rLv,rUc=rUc,rLc=rLc,rUcv=rUcv,rLcv=rLcv,epsilon=epsilon,alpha=0.01,p=p,omega=omega,omegav=omegav,mu=mu,f=f)
        out0 <- ode(states, times=c(0,0), func = "derivs", parms = parameters0,dllname = "SVIC", initfunc = "initmod",nout=32,outnames=paste0("out",0:31),method = "lsoda")
        get.matricies(out0)
        Re.alpha.delta.start<-getR0(Fmat,Vmat)
        
        parameters1<-c(b1=b1,b2=b2,gamma=gamma,rUv=rUv,rLv=rLv,rUc=rUc,rLc=rLc,rUcv=rUcv,rLcv=rLcv,epsilon=epsilon,alpha=alpha1,p=p,omega=omega,omegav=omegav,mu=mu,f=f)
        out1 <- ode(states, times=times, func = "derivs", parms = parameters1,dllname = "SVIC", initfunc = "initmod",nout=32,outnames=paste0("out",0:31),method = "radau")
        epi.equi.states<-out1[nrow(out1),c("S","V","I_0","I_V","I_C","I_C_V","C","C_V")]
        
        out2 <- ode(epi.equi.states, times=c(0,0), func = "derivs", parms = parameters0,dllname = "SVIC", initfunc = "initmod",nout=32,outnames=paste0("out",0:31),method = "lsoda")
        get.matricies(out2)
        Re.alpha.delta.epi.equi<-getR0(Fmat,Vmat)
        
        for (alpha2 in virulence.steps.mod)
        {          
          parameters2<-c(b1=b1,b2=b2,gamma=gamma,rUv=rUv,rLv=rLv,rUc=rUc,rLc=rLc,rUcv=rUcv,rLcv=rLcv,epsilon=epsilon,alpha=alpha2,p=p,omega=omega,omegav=omegav,mu=mu,f=f)
          out3 <- ode(epi.equi.states, times=c(0,0), func = "derivs", parms = parameters2,dllname = "SVIC", initfunc = "initmod",nout=32,outnames=paste0("out",0:31))
          get.matricies(out3)
          RE.invader<-getR0(Fmat,Vmat)
          RE.invader.vec<-c(RE.invader.vec,RE.invader)
        }
        print(paste0("finished alpha1 = ",alpha1))
      }
      RE.invader.mat<-matrix(RE.invader.vec,length(virulence.steps.mod),length(virulence.steps.mod),byrow = T)
      colnames(RE.invader.mat)<-virulence.steps.mod
      rownames(RE.invader.mat)<-virulence.steps.mod
      image(RE.invader.mat>=1)
      ess.result<-ess.analysis(RE.invader.mat)
      data[[index]]<-c(
                     "rUv"=rUv,
                     "rLv"=rLv,
                     "Re.alpha.delta.start"=Re.alpha.delta.start,
                     "Re.alpha.delta.epi.equi"=Re.alpha.delta.epi.equi,
                     "pip.motif"=na.omit(ess.result[c(1,3,4)]),
                     "alpha.ess"=ess.result[2]
                   )
      print(paste0("finished index = ",index))
      index<-index+1
    }
  }
  saveRDS(data,file="~/Documents/GitHub/covid_vaccine_evo/sim.data/rUc0.5rLc0.75p.vacc0.1alpha.optim0.00875.RDS")
}

#plot.mat.R0.obs<-matrix(R0.obs.vec,res,res,byrow = T) #populate matricies
#plot.mat.R0.mutant<-matrix(R0.mutant.vec,res,res,byrow = T) #populate matricies
#.mat<-plot.mat.R0.mutant-plot.mat.R0.obs
#plot.s(s.mat,s.colors[[color.index]],s.col.vals[[color.index]])
#contour(plot.mat.R0.mutant-plot.mat.R0.obs,add=T)
#mtext(expression('r'[L]),side = 2,line=2.5)
#mtext("10% vaccinated",line=2,cex=1.25)

# 50% vacc
{
  get.states(p.C=.25,p.I=0,p.vacc=.5)
  plot.mat.R0.obs<-matrix(NA,res,res) #build matricies to populate
  plot.mat.R0.mutant<-matrix(NA,res,res) #build matricies to populate
  R0.obs.vec<-c()
  R0.mutant.vec<-c()
  
  for (rUx in seq(0,1,length.out = res))
  {
    for (rLx in seq(0,1,length.out = res))
    {
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=vir.obs,prop=prop)
      out2 <- ode(states, times=c(0,0), func = "derivs", parms = parameters,
                  dllname = "SVIC", initfunc = "initmod",nout=18,outnames=paste0("out",0:17))
      get.matricies(out2)
      R0.obs<-getR0(Fmat,Vmat)
      
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=1.5*vir.obs,prop=prop)
      out2 <- ode(states, times=c(0,0), func = "derivs", parms = parameters,
                  dllname = "SVIC", initfunc = "initmod",nout=18,outnames=paste0("out",0:17))
      get.matricies(out2)
      R0.mutant<-getR0(Fmat,Vmat)
      
      R0.obs.vec<-c(R0.obs.vec,R0.obs)
      R0.mutant.vec<-c(R0.mutant.vec,R0.mutant)
    }
  }
  
  plot.mat.R0.obs<-matrix(R0.obs.vec,res,res,byrow = T) #populate matricies
  plot.mat.R0.mutant<-matrix(R0.mutant.vec,res,res,byrow = T) #populate matricies
  s.mat<-plot.mat.R0.mutant-plot.mat.R0.obs
  plot.s(s.mat,s.colors[[color.index]],s.col.vals[[color.index]])
  contour(plot.mat.R0.mutant-plot.mat.R0.obs,add=T)
  mtext("50% vaccinated",line=2,cex=1.25)
  #mtext(expression('selection for '*alpha*' = 0.01'),side=3,line=4,font=2,cex=1.2)
}

# 90% vacc
{
  get.states(p.C=.25,p.I=0,p.vacc=.9)
  plot.mat.R0.obs<-matrix(NA,res,res) #build matricies to populate
  plot.mat.R0.mutant<-matrix(NA,res,res) #build matricies to populate
  R0.obs.vec<-c()
  R0.mutant.vec<-c()
  
  for (rUx in seq(0,1,length.out = res))
  {
    for (rLx in seq(0,1,length.out = res))
    {
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=vir.obs,prop=prop)
      out2 <- ode(states, times=c(0,0), func = "derivs", parms = parameters,
                  dllname = "SVIC", initfunc = "initmod",nout=18,outnames=paste0("out",0:17))
      get.matricies(out2)
      R0.obs<-getR0(Fmat,Vmat)
      
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=1.5*vir.obs,prop=prop)
      out2 <- ode(states, times=c(0,0), func = "derivs", parms = parameters,
                  dllname = "SVIC", initfunc = "initmod",nout=18,outnames=paste0("out",0:17))
      get.matricies(out2)
      R0.mutant<-getR0(Fmat,Vmat)
      
      R0.obs.vec<-c(R0.obs.vec,R0.obs)
      R0.mutant.vec<-c(R0.mutant.vec,R0.mutant)
    }
  }
  
  plot.mat.R0.obs<-matrix(R0.obs.vec,res,res,byrow = T) #populate matricies
  plot.mat.R0.mutant<-matrix(R0.mutant.vec,res,res,byrow = T) #populate matricies
  s.mat<-plot.mat.R0.mutant-plot.mat.R0.obs
  plot.s(s.mat,s.colors[[color.index]],s.col.vals[[color.index]])
  contour(plot.mat.R0.mutant-plot.mat.R0.obs,add=T)
  mtext("90% vaccinated",line=2,cex=1.25)
}


### optim vir = alpha[B.1.1.7]
par(mar=c(2,2,2,2))
optim.vir.assumed<-.0075 #set to either .005, .01, .0025
color.index<-2 #index to match colors and color values to right analysis. 1 -> optim.vir.assumed=.01, 2 -> optim.vir.assumed=.005, 3 -> optim.vir.assumed=.0025

rU<-.0 #vaccinated class
rL<-0 #vaccinated class
rUc<-0 #convalescent class
rLc<-0 #convalescent class

get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,200),vr=vir.obs,b2=b2,tol=1e-10)$root

rUc<-rUc.fix
rLc<-rLc.fix

# 10% vacc
{
  get.states(p.C=.25,p.I=0,p.vacc=.1)
  plot.mat.R0.obs<-matrix(NA,res,res) #build matricies to populate
  plot.mat.R0.mutant<-matrix(NA,res,res) #build matricies to populate
  R0.obs.vec<-c()
  R0.mutant.vec<-c()
  for (rUx in seq(0,1,length.out = res))
  {
    for (rLx in seq(0,1,length.out = res))
    {
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=vir.obs,prop=prop)
      out2 <- ode(states, times=c(0,0), func = "derivs", parms = parameters,
                  dllname = "SVIC", initfunc = "initmod",nout=18,outnames=paste0("out",0:17))
      get.matricies(out2)
      R0.obs<-getR0(Fmat,Vmat)
      
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=1.5*vir.obs,prop=prop)
      out2 <- ode(states, times=c(0,0), func = "derivs", parms = parameters,
                  dllname = "SVIC", initfunc = "initmod",nout=18,outnames=paste0("out",0:17))
      get.matricies(out2)
      R0.mutant<-getR0(Fmat,Vmat)
      
      R0.obs.vec<-c(R0.obs.vec,R0.obs)
      R0.mutant.vec<-c(R0.mutant.vec,R0.mutant)
    }
  }
  
  plot.mat.R0.obs<-matrix(R0.obs.vec,res,res,byrow = T) #populate matricies
  plot.mat.R0.mutant<-matrix(R0.mutant.vec,res,res,byrow = T) #populate matricies
  s.mat<-plot.mat.R0.mutant-plot.mat.R0.obs
  plot.s(s.mat,s.colors[[color.index]],s.col.vals[[color.index]])
  contour(plot.mat.R0.mutant-plot.mat.R0.obs,add=T)
  #mtext(expression('r'[L]),side = 2,line=2.5)
  mtext(expression(alpha['optim']*' = '*alpha[B.1.1.7]),side=2,line=7,cex=1.25)
}

# 50% vacc
{
  get.states(p.C=.25,p.I=0,p.vacc=.5)
  plot.mat.R0.obs<-matrix(NA,res,res) #build matricies to populate
  plot.mat.R0.mutant<-matrix(NA,res,res) #build matricies to populate
  R0.obs.vec<-c()
  R0.mutant.vec<-c()
  
  for (rUx in seq(0,1,length.out = res))
  {
    for (rLx in seq(0,1,length.out = res))
    {
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=vir.obs,prop=prop)
      out2 <- ode(states, times=c(0,0), func = "derivs", parms = parameters,
                  dllname = "SVIC", initfunc = "initmod",nout=18,outnames=paste0("out",0:17))
      get.matricies(out2)
      R0.obs<-getR0(Fmat,Vmat)
      
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=1.5*vir.obs,prop=prop)
      out2 <- ode(states, times=c(0,0), func = "derivs", parms = parameters,
                  dllname = "SVIC", initfunc = "initmod",nout=18,outnames=paste0("out",0:17))
      get.matricies(out2)
      R0.mutant<-getR0(Fmat,Vmat)
      
      R0.obs.vec<-c(R0.obs.vec,R0.obs)
      R0.mutant.vec<-c(R0.mutant.vec,R0.mutant)
    }
  }
  
  plot.mat.R0.obs<-matrix(R0.obs.vec,res,res,byrow = T) #populate matricies
  plot.mat.R0.mutant<-matrix(R0.mutant.vec,res,res,byrow = T) #populate matricies
  s.mat<-plot.mat.R0.mutant-plot.mat.R0.obs
  plot.s(s.mat,s.colors[[color.index]],s.col.vals[[color.index]])
  contour(plot.mat.R0.mutant-plot.mat.R0.obs,add=T)
  #mtext(expression('selection for '*alpha*' = 0.01'),side=3,line=4,font=2,cex=1.2)
}

# 90% vacc
{
  get.states(p.C=.25,p.I=0,p.vacc=.9)
  plot.mat.R0.obs<-matrix(NA,res,res) #build matricies to populate
  plot.mat.R0.mutant<-matrix(NA,res,res) #build matricies to populate
  R0.obs.vec<-c()
  R0.mutant.vec<-c()
  
  for (rUx in seq(0,1,length.out = res))
  {
    for (rLx in seq(0,1,length.out = res))
    {
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=vir.obs,prop=prop)
      out2 <- ode(states, times=c(0,0), func = "derivs", parms = parameters,
                  dllname = "SVIC", initfunc = "initmod",nout=18,outnames=paste0("out",0:17))
      get.matricies(out2)
      R0.obs<-getR0(Fmat,Vmat)
      
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=1.5*vir.obs,prop=prop)
      out2 <- ode(states, times=c(0,0), func = "derivs", parms = parameters,
                  dllname = "SVIC", initfunc = "initmod",nout=18,outnames=paste0("out",0:17))
      get.matricies(out2)
      R0.mutant<-getR0(Fmat,Vmat)
      
      R0.obs.vec<-c(R0.obs.vec,R0.obs)
      R0.mutant.vec<-c(R0.mutant.vec,R0.mutant)
    }
  }
  
  plot.mat.R0.obs<-matrix(R0.obs.vec,res,res,byrow = T) #populate matricies
  plot.mat.R0.mutant<-matrix(R0.mutant.vec,res,res,byrow = T) #populate matricies
  s.mat<-plot.mat.R0.mutant-plot.mat.R0.obs
  plot.s(s.mat,s.colors[[color.index]],s.col.vals[[color.index]])
  contour(plot.mat.R0.mutant-plot.mat.R0.obs,add=T)
}


### optim vir = alpha[ansc] * 1.25
par(mar=c(2,2,2,2))
optim.vir.assumed<-.00625 #set to either .005, .01, .0025
color.index<-3 #index to match colors and color values to right analysis. 1 -> optim.vir.assumed=.01, 2 -> optim.vir.assumed=.005, 3 -> optim.vir.assumed=.0025

rU<-.0 #vaccinated class
rL<-0 #vaccinated class
rUc<-0 #convalescent class
rLc<-0 #convalescent class

get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,200),vr=vir.obs,b2=b2,tol=1e-10)$root

rUc<-rUc.fix
rLc<-rLc.fix

# 10% vacc
{
  get.states(p.C=.25,p.I=0,p.vacc=.1)
  plot.mat.R0.obs<-matrix(NA,res,res) #build matricies to populate
  plot.mat.R0.mutant<-matrix(NA,res,res) #build matricies to populate
  R0.obs.vec<-c()
  R0.mutant.vec<-c()
  
  for (rUx in seq(0,1,length.out = res))
  {
    for (rLx in seq(0,1,length.out = res))
    {
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=vir.obs,prop=prop)
      out2 <- ode(states, times=c(0,0), func = "derivs", parms = parameters,
                  dllname = "SVIC", initfunc = "initmod",nout=18,outnames=paste0("out",0:17))
      get.matricies(out2)
      R0.obs<-getR0(Fmat,Vmat)
      
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=1.5*vir.obs,prop=prop)
      out2 <- ode(states, times=c(0,0), func = "derivs", parms = parameters,
                  dllname = "SVIC", initfunc = "initmod",nout=18,outnames=paste0("out",0:17))
      get.matricies(out2)
      R0.mutant<-getR0(Fmat,Vmat)
      
      R0.obs.vec<-c(R0.obs.vec,R0.obs)
      R0.mutant.vec<-c(R0.mutant.vec,R0.mutant)
    }
  }
  
  plot.mat.R0.obs<-matrix(R0.obs.vec,res,res,byrow = T) #populate matricies
  plot.mat.R0.mutant<-matrix(R0.mutant.vec,res,res,byrow = T) #populate matricies
  s.mat<-plot.mat.R0.mutant-plot.mat.R0.obs
  plot.s(s.mat,s.colors[[color.index]],s.col.vals[[color.index]])
  contour(plot.mat.R0.mutant-plot.mat.R0.obs,add=T)
  mtext(expression(alpha['optim']*' = 1.25*'*alpha[ansc]),side=2,line=7,cex=1.25)
}

# 50% vacc
{
  get.states(p.C=.25,p.I=0,p.vacc=.5)
  plot.mat.R0.obs<-matrix(NA,res,res) #build matricies to populate
  plot.mat.R0.mutant<-matrix(NA,res,res) #build matricies to populate
  R0.obs.vec<-c()
  R0.mutant.vec<-c()
  
  for (rUx in seq(0,1,length.out = res))
  {
    for (rLx in seq(0,1,length.out = res))
    {
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=vir.obs,prop=prop)
      out2 <- ode(states, times=c(0,0), func = "derivs", parms = parameters,
                  dllname = "SVIC", initfunc = "initmod",nout=18,outnames=paste0("out",0:17))
      get.matricies(out2)
      R0.obs<-getR0(Fmat,Vmat)
      
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=1.5*vir.obs,prop=prop)
      out2 <- ode(states, times=c(0,0), func = "derivs", parms = parameters,
                  dllname = "SVIC", initfunc = "initmod",nout=18,outnames=paste0("out",0:17))
      get.matricies(out2)
      R0.mutant<-getR0(Fmat,Vmat)
      
      R0.obs.vec<-c(R0.obs.vec,R0.obs)
      R0.mutant.vec<-c(R0.mutant.vec,R0.mutant)
    }
  }
  
  plot.mat.R0.obs<-matrix(R0.obs.vec,res,res,byrow = T) #populate matricies
  plot.mat.R0.mutant<-matrix(R0.mutant.vec,res,res,byrow = T) #populate matricies
  s.mat<-plot.mat.R0.mutant-plot.mat.R0.obs
  plot.s(s.mat,s.colors[[color.index]],s.col.vals[[color.index]])
  contour(plot.mat.R0.mutant-plot.mat.R0.obs,add=T)
}

# 90% vacc
{
  get.states(p.C=.25,p.I=0,p.vacc=.9)
  plot.mat.R0.obs<-matrix(NA,res,res) #build matricies to populate
  plot.mat.R0.mutant<-matrix(NA,res,res) #build matricies to populate
  R0.obs.vec<-c()
  R0.mutant.vec<-c()
  
  for (rUx in seq(0,1,length.out = res))
  {
    for (rLx in seq(0,1,length.out = res))
    {
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=vir.obs,prop=prop)
      out2 <- ode(states, times=c(0,0), func = "derivs", parms = parameters,
                  dllname = "SVIC", initfunc = "initmod",nout=18,outnames=paste0("out",0:17))
      get.matricies(out2)
      R0.obs<-getR0(Fmat,Vmat)
      
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=1.5*vir.obs,prop=prop)
      out2 <- ode(states, times=c(0,0), func = "derivs", parms = parameters,
                  dllname = "SVIC", initfunc = "initmod",nout=18,outnames=paste0("out",0:17))
      get.matricies(out2)
      R0.mutant<-getR0(Fmat,Vmat)
      
      R0.obs.vec<-c(R0.obs.vec,R0.obs)
      R0.mutant.vec<-c(R0.mutant.vec,R0.mutant)
    }
  }
  
  plot.mat.R0.obs<-matrix(R0.obs.vec,res,res,byrow = T) #populate matricies
  plot.mat.R0.mutant<-matrix(R0.mutant.vec,res,res,byrow = T) #populate matricies
  s.mat<-plot.mat.R0.mutant-plot.mat.R0.obs
  plot.s(s.mat,s.colors[[color.index]],s.col.vals[[color.index]])
  contour(plot.mat.R0.mutant-plot.mat.R0.obs,add=T)
}

#legend
{
  par(mar=c(2,1,2,5))
  yy<-seq(0,length(s.col.vals[[color.index]]),1)
  plot(0,0,type="n",xlim=c(0,1),ylim=c(.5,length(s.col.vals[[color.index]])+.5),xlab="",ylab="",axes=F)
  color.legend(0,0,1,length(s.col.vals[[color.index]]),legend=NULL,s.colors[[color.index]],gradient="y")
  axis(4,at=seq(0,length(s.col.vals[[color.index]])-1,length.out = 5)+1,labels = s.col.vals[[color.index]][seq(0,length(s.col.vals[[color.index]])-1,length.out = 5)+1])
  mtext("selection",side=4,line =2)
}

mtext(expression('lower respiratory tract protection (r'["L,V"]*')'),side = 2,line=1.5,cex=1.5,outer=T)
mtext(expression('upper respiratory tract protection (r'["U,V"]*')'),side = 1,line=2,cex=1.5,outer=T,adj=3/7)

## copy to clipboard with width of 1051, height of 843


