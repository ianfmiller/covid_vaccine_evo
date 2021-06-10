#### set up
source("covid_SVIC_functions.R")

library(deSolve)
library(rootSolve)
library(plotrix)

system("R CMD SHLIB SVIC.c")
dyn.load(paste("SVIC", .Platform$dynlib.ext, sep = ""))

#### define model parameters

gamma<-1

optim.vir.assumed<-.0075 #set to either .01, .005, .0025

vir.obs<-.0075

prop<-66.666

R0.assumed<-3.75

rU<-.0 #vaccinated class
rL<-0 #vaccinated class
rUc<-0 #convalescent class
rLc<-0 #convalescent class

frac_lower<-.5 # % contribution of lower respiratory infection to overall transmission


### set trade-off (b1 + b2) according to assumed ES virulence, observed virulence, and R0 at observed virulence

## set trade-off shape (b2 only) according to assumed ES virulence

get.states(0,0,0)

b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root

## set trade-off scaling (b1) so that R0=2.5 at vir.obs

b1<-uniroot(R0.search,c(0,100),vr=vir.obs,b2=b2,tol=1e-10)$root

## checks for setting b1 and b2

#check to ensure that fit b2 is effectively independent of b1
#uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-10)$root
#uniroot(b2.search,c(0,1),b1=2,optim.vir.assumed=optim.vir.assumed,tol=1e-10)$root
#uniroot(b2.search,c(0,1),b1=3,optim.vir.assumed=optim.vir.assumed,tol=1e-10)$root

#check to ensure that b1 and b2 give desired R0 at vir.obs
#find.R0(vir.obs,b1,b2)

#check that optim.vir.assumed is in fact optimal given b1 and b2, and that b1 and b2 and give desired R0 at vir.obs

#vs<-seq(0,.02,.00001)
#rs<-unlist(lapply(vs,find.R0,b1=b1,b2=b2))
#plot(vs,rs,ylim=c(3.74,3.76),xlab=expression(alpha),ylab=expression(R[0]))
#abline(v=optim.vir.assumed)
#abline(h=max(rs))
#abline(v=vir.obs)
#abline(h=R0.assumed)


#### calculate selection coefficient for combinations of rL and rU

### plotting tools

col.trans<-.6

s.blues1<-hsv(.666,1,1,seq(1,.001,length.out = 400)^col.trans)
s.reds1<-hsv(1,1,1,seq(.001,1,length.out = 400)^col.trans)
s.colors1<-c(s.blues1,"white",s.reds1)
s.col.vals1<-seq(-.55,.55,length.out = 801)

s.blues2<-hsv(.666,1,1,seq(1,.001,length.out = 400)^col.trans)
s.reds2<-hsv(1,1,1,seq(.001,1,length.out = 400)^col.trans)
s.colors2<-c(s.blues2,"white",s.reds2)
s.col.vals2<-seq(-.55,.55,length.out = 801)

s.blues3<-hsv(.666,1,1,seq(1,.001,length.out = 400)^col.trans)
s.reds3<-hsv(1,1,1,seq(.001,1,length.out = 400)^col.trans)
s.colors3<-c(s.blues3,"white",s.reds3)
s.col.vals3<-seq(-.55,.55,length.out = 801)

s.colors<-list(s.colors1,s.colors2,s.colors3)
s.col.vals<-list(s.col.vals1,s.col.vals2,s.col.vals3)

res<-101

### set plot window

layout(matrix(c(1,1,2,2,3,3,10,10,4,4,5,5,6,6,10,10,7,7,8,8,9,9,10,10),3,8,byrow = T))
par(mar=c(2,2,2,2),oma=c(4,8,4,0))

#### analysis

rUc.fix<-.25 #convalescent class
rLc.fix<-.5 #convalescent class

### optim vir = 2 * obs vir

optim.vir.assumed<-.0075*2 #set to either .005, .01, .0025
color.index<-1 #index to match colors and color values to right analysis. 1 -> optim.vir.assumed=.01, 2 -> optim.vir.assumed=.005, 3 -> optim.vir.assumed=.0025

rU<-.0 #vaccinated class
rL<-0 #vaccinated class
rUc<-0 #convalescent class
rLc<-0 #convalescent class

get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,100),vr=vir.obs,b2=b2,tol=1e-10)$root

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
      
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=2*vir.obs,prop=prop)
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
  plot.outcome(plot.mat.R0.obs,plot.mat.R0.mutant)
  #mtext(expression('r'[L]),side = 2,line=2.5)
  mtext("10% vaccinated",line=2,cex=1.25)
  mtext(expression(alpha['optim']*' = 2*'*alpha[B.1.1.7]),side=2,line=7,cex=1.25)
}

# 20% vacc
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
      
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=2*vir.obs,prop=prop)
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
  plot.outcome(plot.mat.R0.obs,plot.mat.R0.mutant)
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
      
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=2*vir.obs,prop=prop)
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
  plot.outcome(plot.mat.R0.obs,plot.mat.R0.mutant)
  mtext("90% vaccinated",line=2,cex=1.25)
}

### optim vir = 1.5 * obs vir
par(mar=c(2,2,2,2))
optim.vir.assumed<-.0075 #set to either .005, .01, .0025
color.index<-2 #index to match colors and color values to right analysis. 1 -> optim.vir.assumed=.01, 2 -> optim.vir.assumed=.005, 3 -> optim.vir.assumed=.0025

rU<-.0 #vaccinated class
rL<-0 #vaccinated class
rUc<-0 #convalescent class
rLc<-0 #convalescent class

get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,100),vr=vir.obs,b2=b2,tol=1e-10)$root

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
      
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=2*vir.obs,prop=prop)
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
  plot.outcome(plot.mat.R0.obs,plot.mat.R0.mutant)
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
      
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=2*vir.obs,prop=prop)
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
  plot.outcome(plot.mat.R0.obs,plot.mat.R0.mutant)
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
      
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=2*vir.obs,prop=prop)
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
  plot.outcome(plot.mat.R0.obs,plot.mat.R0.mutant)
}

### optim vir = obs vir
par(mar=c(2,2,2,2))
optim.vir.assumed<-.00625 #set to either .005, .01, .0025
color.index<-3 #index to match colors and color values to right analysis. 1 -> optim.vir.assumed=.01, 2 -> optim.vir.assumed=.005, 3 -> optim.vir.assumed=.0025

rU<-.0 #vaccinated class
rL<-0 #vaccinated class
rUc<-0 #convalescent class
rLc<-0 #convalescent class

get.states(0,0,0)
b2<-uniroot(b2.search,c(0,1),b1=1,optim.vir.assumed=optim.vir.assumed,tol=1e-15)$root
b1<-uniroot(R0.search,c(0,100),vr=vir.obs,b2=b2,tol=1e-10)$root

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
      
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=2*vir.obs,prop=prop)
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
  plot.outcome(plot.mat.R0.obs,plot.mat.R0.mutant)
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
      
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=2*vir.obs,prop=prop)
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
  plot.outcome(plot.mat.R0.obs,plot.mat.R0.mutant)
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
      
      parameters<-c(b1=b1,b2=b2,gamma=gamma,rU=rUx,rL=rLx,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=2*vir.obs,prop=prop)
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
  plot.outcome(plot.mat.R0.obs,plot.mat.R0.mutant)
}

mtext(expression('lower respiratory tract protection (r'["L,V"]*')'),side = 2,line=1.5,cex=1.5,outer=T)
mtext(expression('upper respiratory tract protection (r'["U,V"]*')'),side = 1,line=2,cex=1.5,outer=T,adj=3/7)

#legend
plot(0,0,type="n",xlab="",ylab="",bty="n",axes = F,xlim=c(0,1),ylim=c(0,1))
par(oma=c(4,0,4,0),mar=c(5,0,4,0))
legend("left",yjust=.5,
       legend=c(
         expression(' 1 > '*R[E](alpha[B.1.1.7])*' > '*R[E](1.5*'*'*alpha[B.1.1.7])),NA,NA,
         expression(' 1 > '*R[E](1.5*'*'*alpha[B.1.1.7])*' > '*R[E](alpha[B.1.1.7])),NA,NA,
         expression(' '*R[E](alpha[B.1.1.7])*' > 1 & '),expression(R[E](alpha[B.1.1.7])*' > '*R[E](1.5*'*'*alpha[B.1.1.7])),NA,
         expression(' '*R[E](1.5*'*'*alpha[B.1.1.7])*' > 1 > '*R[E](alpha[B.1.1.7])),NA,NA,
         expression(' '*R[E](1.5*'*'*alpha[B.1.1.7])*' > 1 &'),expression(R[E](1.5*'*'*alpha[B.1.1.7])*' > '*R[E](alpha[B.1.1.7])),NA
       ),
       col=c(
         viridis(4)[1],NA,NA,
         viridis(4,alpha=.5)[1],NA,NA,
         viridis(4)[2],NA,NA,
         viridis(4)[3],NA,NA,
         viridis(4)[4],NA,NA
       ),
       pch=15,cex=1.1,pt.cex = 4,bty="n",y.intersp = 1.5)

## copy to clipboard with width of 1185, height of 843