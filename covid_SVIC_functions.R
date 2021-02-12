### convenience functions for analysis

library(viridis)

## functions for setting b1 according to assumed ES virulence, dr, and gamma

get.matricies<-function(output) #calculates F and V matricies (used in next-gen R0 calculation) from output of SIR model
{
  Fmat<<-matrix(output[1,8:16],3,3,byrow = T)
  Vmat<<-matrix(output[1,17:25],3,3,byrow = T)
}

getR0<-function(Fmat,Vmat) #calculates R0 from F and V matricies
{
  values<-eigen(Fmat %*% solve(Vmat))$values
  R0<<-0
  for(i in 1:length(values))
  {
    if(Im(values)[i]==0) {R0<<-max(R0,Re(values)[i])}
  }
  R0
}

get.states<-function(p.C, p.I, p.vacc) #set initial conditions
{
  S_0=1-p.vacc-p.I-p.C #susceptible
  V_0=p.vacc #vaccinated
  I_n_0=p.I #infected--naive--resdient strain
  I_v_0=0 #infected--vaccinated--resdient strain
  I_c_0=0 #infected--convalescent--resident strain
  C_0=p.C #convalescent
  
  states<<-c(S=S_0,
            V=V_0,
            I_n=I_n_0,
            I_v=I_v_0,
            I_c=I_c_0,
            C=C_0)
}


find.R0<-function(v,b1,b2) # get R0 given vr, b1, b2
{
  parameters <- c(b1=b1,b2=b2,gamma=gamma,rU=rU,rL=rL,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=v,prop=prop)
  new.out <- ode(states, c(0,0), func = "derivs", parms = parameters,
                 dllname = "SVIC", initfunc = "initmod",nout=18,outnames=paste0("out",0:17))
  get.matricies(new.out)
  R0.calc<-getR0(Fmat,Vmat)
  R0.calc
}

R0.search<-function(vr,b1,b2) # get differen between assumed R0 and R0 given vr, b1, b2
{
  R0.assumed-find.R0(vr,b1,b2)
}

b2.search<-function(b1,b2,optim.vir.assumed) #get abs difference between optim vir assumed and true optim vir given b1,b2
{
  optim.vir.assumed-optimize(find.R0,b1=b1,b2=b2,interval = c(0.0025,1),maximum = T,tol=1e-10)$maximum
}

find.optim.vir<-function(vsteps,b1,b2,rU,rL,rUc,rLc) # get optim vir given b1, b2
{
  rU<-rU
  rL<-rL
  rUc<-rUc
  rLc<-rLc
  
  R0s<-c()
  for(v in vsteps)
  {
    parameters <- c(b1=b1,b2=b2,gamma=gamma,rU=rU,rL=rL,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=v,prop=prop)
    new.out <- ode(states, c(0,0), func = "derivs", parms = parameters,
                   dllname = "SVIC", initfunc = "initmod",nout=18,outnames=paste0("out",0:17))
    get.matricies(new.out)
    R0.calc<-getR0(Fmat,Vmat)
    R0s<-c(R0s,R0.calc)
  }
  list(vsteps[which.max(R0s)],R0s)
}

plot.s<-function(plot.mat,cols,col.vals) #plotting function
{
  plot(0,0,type="n",xlim=c(-(1/(res-1))/2,1+(1/(res-1))/2),ylim=c(-(1/(res-1))/2,1+(1/(res-1))/2),xlab=expression('r'[U]),ylab=expression('r'[L]),cex.lab=2)
  xx<-seq(0,1,length.out = res)
  yy<-seq(0,1,length.out = res)
  for(i in 1:res)
  {
    for(j in 1:res)
    {
      rect(xx[i]-(1/(res-1))/2,yy[j]-(1/(res-1))/2,xx[i]+(1/(res-1))/2,yy[j]+(1/(res-1))/2,col = cols[which.min(abs(plot.mat[i,j]-col.vals))],border=NA)
    }
  }
}

outcome.col.func<-function(R0.obs,R0.mutant)
{
  if(abs(R0.obs)<1e-10) {R0.obs<-0}
  if(abs(R0.mutant)<1e-10) {R0.mutant<-0}
  if(R0.obs<R0.mutant) # selection for increased virulence
  {
    if(R0.obs<1 && R0.mutant<1) {val<-viridis(4,alpha=.75)[1]} # erad w/ selection for increased virulence
    if(R0.obs<1 && R0.mutant>=1) {val<-viridis(4)[3]} # erad w/ evol escape
    if(R0.obs>=1 && R0.mutant>=1) {val<-viridis(4)[4]} # selection for increased virulence
  }
  
  if(R0.obs>=R0.mutant) # selection against increased virulence
  {
    if(R0.obs<1) {val<-viridis(4)[1]} # erad w/ selection against increased virulence
    if(R0.obs>=1) {val<-viridis(4)[2]} # selection against increased virulence, no erad
  }
  return(val)
}

plot.outcome<-function(plot.mat.R0.obs,plot.mat.R0.mutatnt)
{
  plot(0,0,type="n",xlim=c(-(1/(res-1))/2,1+(1/(res-1))/2),ylim=c(-(1/(res-1))/2,1+(1/(res-1))/2),xlab=expression('r'[U]),ylab=expression('r'[L]),cex.lab=2)
  xx<-seq(0,1,length.out = res)
  yy<-seq(0,1,length.out = res)
  for(i in 1:res)
  {
    for(j in 1:res)
    {
      rect(xx[i]-(1/(res-1))/2,yy[j]-(1/(res-1))/2,xx[i]+(1/(res-1))/2,yy[j]+(1/(res-1))/2,col = outcome.col.func(plot.mat.R0.obs[i,j],plot.mat.R0.mutant[i,j]),border=NA)
    }
  }
}






