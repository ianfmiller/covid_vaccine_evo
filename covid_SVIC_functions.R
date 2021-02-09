### convenience functions for analysis

## functions for setting b1 according to assumed ES virulence, dr, and gamma

get.matricies<-function(output) #calculates F and V matricies (used in next-gen R0 calculation) from output of SIR model
{
  Fmat<<-matrix(output[1,11:46],6,6,byrow = T)
  Vmat<<-matrix(output[1,47:82],6,6,byrow = T)
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
  I_pn_0=p.I #infected--presymptomatic--naive
  I_pv_0=0 #infected--presymptomatic--vaccinated
  I_pc_0=0 #infected--presymptomatic--convalescent
  I_sn_0=p.I #infected--symptomatic--naive
  I_sv_0=0 #infected--symptomatic--vaccinated
  I_sc_0=0 #infected--symptomatic--convalescent
  C_0=p.C #convalescent
  
  states<<-c(S=S_0,
            V=V_0,
            I_pn=I_pn_0,
            I_pv=I_pv_0,
            I_pc=I_pc_0,
            I_sn=I_sn_0,
            I_sv=I_sv_0,
            I_sc=I_sc_0,
            C=C_0)
}


find.R0<-function(v,b1,b2) # get R0 given vr, b1, b2
{
  parameters <- c(b1=b1,b2=b2,gamma=gamma,sigma=sigma,rUv=rUv,rLv=rLv,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=v,x=x)
  new.out <- ode(states, c(0,0), func = "derivs", parms = parameters,
                 dllname = "SVIC", initfunc = "initmod",nout=72,outnames=paste0("out",0:71))
  get.matricies(new.out)
  R0.calc<-getR0(Fmat,Vmat)
  R0.calc
}

R0.search<-function(vr,b1,b2) # get difference between assumed R0 and R0 given vr, b1, b2
{
  R0.assumed-find.R0(vr,b1,b2)
}

b2.search<-function(b1,b2,optim.vir.assumed) #get abs difference between optim vir assumed and true optim vir given b1,b2
{
  optim.vir.assumed-optimize(find.R0,b1=b1,b2=b2,interval = c(0,.005),maximum = T,tol=1e-15)$maximum
}

find.optim.vir<-function(vsteps,b1,b2,rUv,rLv,rUc,rLc) # get optim vir given b1, b2
{
  rUv<-rUv
  rLv<-rLv
  rUc<-rUc
  rLc<-rLc
  
  R0s<-c()
  for(v in vsteps)
  {
    parameters <- c(b1=b1,b2=b2,gamma=gamma,sigma=sigma,rUv=rUv,rLv=rLv,rUc=rUc,rLc=rLc,frac_lower=frac_lower,v=v,x=x)
    new.out <- ode(states, c(0,0), func = "derivs", parms = parameters,
                   dllname = "SVIC", initfunc = "initmod",nout=72,outnames=paste0("out",0:71))
    get.matricies(new.out)
    R0.calc<-getR0(Fmat,Vmat)
    R0s<-c(R0s,R0.calc)
  }
  list(vsteps[which.max(R0s)],R0s)
}

plot.s<-function(plot.mat,cols,col.vals) #plotting function
{
  plot(0,0,type="n",xlim=c(-.025,1.025),ylim=c(-0.025,1.025),xlab=expression('r'[U]),ylab=expression('r'[L]),cex.lab=2)
  xx<-seq(0,1,.05)
  yy<-seq(0,1,.05)
  for(i in 1:21)
  {
    for(j in 1:21)
    {
      rect(xx[i]-.025,yy[j]-.025,xx[i]+.025,yy[j]+.025,col = cols[which.min(abs(plot.mat[i,j]-col.vals))],border=NA)
    }
  }
}
