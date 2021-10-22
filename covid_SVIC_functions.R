### convenience functions for analysis

library(viridis)

## functions for setting b1 according to assumed ES virulence, dr, and gamma

get.matricies<-function(output) #calculates F and V matricies (used in next-gen R0 calculation) from output of SIR model
{
  Fmat<<-matrix(output[1,10:25],4,4,byrow = T)
  Vmat<<-matrix(output[1,26:41],4,4,byrow = T)
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

get.states<-function(p.C, p.C_V, p.I, p.V) #set initial conditions
{
  S_0=1-p.V-p.I-p.C-p.C_V #susceptible
  V_0=p.V #vaccinated
  I_0_0=p.I*(1-p.C-p.C_V-p.V) #infected--naive
  I_V_0=p.I*(p.V) #infected--vaccinated
  I_C_0=p.I*(p.C) #infected--convalescent
  I_C_V_0=p.I*(p.C_V) #infected--convalescent
  C_0=p.C #convalescent
  C_V_0=p.C_V
  
  states<<-c(S=S_0,
            V=V_0,
            I_0=I_0_0,
            I_V=I_V_0,
            I_C=I_C_0,
            I_C_V=I_C_V_0,
            C=C_0,
            C_V=C_V_0)
}


find.R0<-function(alpha,b1,b2) # get R0 given vr, b1, b2
{
  parameters <- c(b1=b1,b2=b2,gamma=gamma,rU=rUv,rL=rLv,rUc=rUc,rLc=rLc,rUcv=rUcv,rLcv=rLcv,epsilon=epsilon,alpha=alpha,p=p,omega=omega,omegav=omegav)
  new.out <- ode(states, c(0,0), func = "derivs", parms = parameters,
                 dllname = "SVIC", initfunc = "initmod",nout=32,outnames=paste0("out",0:31))
  get.matricies(new.out)
  R0.calc<-getR0(Fmat,Vmat)
  R0.calc
}

R0.search<-function(alpha,b1,b2) # get differen between assumed R0 and R0 given vr, b1, b2
{
  R0.assumed-find.R0(alpha,b1,b2)
}

b2.search<-function(b1,b2,optim.alpha.assumed) #get abs difference between optim vir assumed and true optim vir given b1,b2
{
  optim.alpha.assumed-optimize(find.R0,b1=b1,b2=b2,interval = c(0.0025,1),maximum = T,tol=1e-10)$maximum
}

find.optim.vir<-function(vsteps,b1,b2,rU,rL,rUc,rLc) # get optim vir given b1, b2
{
  rU<-rU
  rL<-rL
  rUc<-rUc
  rLc<-rLc
  
  R0s<-c()
  for(alpha in vsteps)
  {
    parameters <- c(b1=b1,b2=b2,gamma=gamma,rU=rUv,rL=rLv,rUc=rUc,rLc=rLc,rUcv=rUcv,rLcv=rLcv,epsilon=epsilon,alpha=alpha,p=p,omega=omega,omegav=omegav)
    new.out <- ode(states, c(0,0), func = "derivs", parms = parameters,
                   dllname = "SVIC", initfunc = "initmod",nout=32,outnames=paste0("out",0:31))
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

plot.simulation<-function(output,legend=T)
{
  S.plot<-output[,"S"]
  V.plot<-output[,"V"]
  I_0.plot<-output[,"I_0"]
  I_V.plot<-output[,"I_V"]
  I_C.plot<-output[,"I_C"]
  I_C_V.plot<-output[,"I_C_V"]
  C.plot<-output[,"C"]
  C_V.plot<-output[,"C_V"]
  plot.x<-1:nrow(output)
  
  plot(plot.x,S.plot,type="l",col="blue",ylim=c(0,1),lwd=2,xlab="days",ylab="frequency")
  points(plot.x,V.plot,type="l",col="darkolivegreen1",lwd=2)
  points(plot.x,I_0.plot,type="l",col="orange",lwd=2)
  points(plot.x,I_V.plot,type="l",col="orangered",lwd=2)
  points(plot.x,I_C.plot,type="l",col="red",lwd=2)
  points(plot.x,I_C_V.plot,type="l",col="violetred",lwd=2)
  points(plot.x,C.plot,type="l",col="darkseagreen1",lwd=2)
  points(plot.x,C_V.plot,type="l",col="darkgreen",lwd=2)
  
  if(legend)
  {
    legend(
      "topright",
      legend=c("S","V","I_0","I_V","I_C","I_C_V","C","C_V"),
      col=c("blue","darkolivegreen1","orange",'orangered',"red","violetred","darkseagreen1","darkgreen"),
      lty=1,
      lwd=2
    )
  }
}

### Scans vector from start to finish to find 1st local maximum
find.peak<-function(x)
{
  peak.index<-NA
  for (z in 4:(length(x)-3)) 
  {
    if (x[z] >= x[z-1] && x[z] >= x[z-2] && x[z] > x[z-3] && x[z] >= x[z+1] && x[z] >= x[z+2] && x[z] > x[z+3]) {peak.index<-z;break}
  }
  return(peak.index)
}

### Used to find virulence strategy with fitness equal to that of local optimum 
find.break.even.point<-function(x,peak.index)
{
  break.even.point<-NA
  for (z in (peak.index+1):(length(x)))
  {
    if (x[z] >= x[peak.index]) {break.even.point<-z; break}
  }
  return(break.even.point)
}

pip.analysis<-function(mat)
{
  output<-rep(NA,length.out=13)
  if(colnames(mat)[1]=="X") {mat<-mat[,-1]}
  mod.mat<-matrix(NA,dim(mat)[1],dim(mat)[2])
  colnames(mod.mat)<-colnames(mat)
  rownames(mod.mat)<-rownames(mat)
  for(j in 1:dim(mod.mat)[1]) {
    for(k in 1:dim(mod.mat)[1]) {
      mod.mat[j,k] <- 1*(mat[j,k]>1)
    }}
  
  row.max<-find.peak(rowSums(mod.mat))
  row.min<-find.peak(-1*rowSums(mod.mat))
  
  col.max<-find.peak(colSums(mod.mat))
  col.min<-find.peak(-1*colSums(mod.mat))
  
  if (is.numeric(col.max)&&is.numeric(row.min)) {if(col.max-row.min<=1) {output[1]<-"ESS"; output[2]<-virulence.steps[col.max]}}
  if (is.numeric(col.min)&&is.numeric(row.max)) {if(col.min-row.max<=1) {output[3]<-"REPELLER"; output[4]<-virulence.steps[row.max]}}
  
  if(isTRUE(as.numeric(output[2])>0)) #used to check for repeller point and get eradication bounds
  {
    up.sub.mat.cols<-which(virulence.steps>max(as.numeric(output[2],output[4])))
    up.sub.mat<-mod.mat[up.sub.mat.cols,up.sub.mat.cols]
    up.sub.mat.col.sums<-colSums(up.sub.mat)
    
    if (length(which(up.sub.mat.col.sums==min(up.sub.mat.col.sums)))>=2)
    {
      upper.bound<-max(which(up.sub.mat.col.sums==min(up.sub.mat.col.sums)))+(dim(mod.mat)[1]-length(up.sub.mat.col.sums))
      lower.bound<-min(which(up.sub.mat.col.sums==min(up.sub.mat.col.sums)))+(dim(mod.mat)[1]-length(up.sub.mat.col.sums))
      if(upper.bound<dim(mat)[1] && is.na(output[3]))
      {
        output[3]<-"REPELLER2"
        output[4]<-virulence.steps[upper.bound+1]
        output[10]<-"mid.erad.lower"
        output[11]<-virulence.steps[lower.bound]
        output[12]<-"mid.erad.upper"
        output[13]<-virulence.steps[upper.bound]
        
      }
      
      if(upper.bound==dim(mat)[1])
      {
        output[6]<-"upper.erad"
        output[7]<-virulence.steps[lower.bound]
      }
    }
    
    low.sub.mat.cols<-which(virulence.steps<output[2])
    low.sub.mat<-mod.mat[low.sub.mat.cols,low.sub.mat.cols]
    low.sub.mat.col.sums<-colSums(low.sub.mat)
    
    if (length(which(low.sub.mat.col.sums==min(low.sub.mat.col.sums)))>=2)
    {
      upper.bound<-max(which(low.sub.mat.col.sums==min(low.sub.mat.col.sums)))
      lower.bound<-min(which(low.sub.mat.col.sums==min(low.sub.mat.col.sums)))
      
      if(lower.bound==1)
      {
        output[8]<-"lower.erad"
        output[9]<-virulence.steps[upper.bound]
      }
    }
  }
  
  
  if (isFALSE(any(diff(colSums(mod.mat))<0)) && any(mod.mat>0)) 
  {
    output[5]<-"selection for hypervirulence"
    if(length(which(colSums(mod.mat)==min(colSums(mod.mat))))>2)
    {
      output[8]<-"lower.erad"
      output[9]<-virulence.steps[max(which(colSums(mod.mat)==min(colSums(mod.mat))))]
    }
  }
  
  
  if (isFALSE(any(diff(colSums(mod.mat))>0)) && any(mod.mat>0)) {output[5]<-"selection for 0 virulence"}
  
  if (!any(is.na(output)==F)) #this loop is for when theres a region around the ESS where nothing can invade. This is due to numerical errors
  {
    intersect(which(colSums(mod.mat)==max(colSums(mod.mat))),which(rowSums(mod.mat)==min(rowSums(mod.mat))))->overlap
    if (length(overlap)>0)
    {
      output[1]<-"ESS"
      output[2]<-mean(virulence.steps[overlap])
      { #get erad bounds
        up.sub.mat.cols<-which(virulence.steps>max(as.numeric(output[2],output[4])))
        up.sub.mat<-mod.mat[up.sub.mat.cols,up.sub.mat.cols]
        up.sub.mat.col.sums<-colSums(up.sub.mat)
        
        if (length(which(up.sub.mat.col.sums==min(up.sub.mat.col.sums)))>=2)
        {
          upper.bound<-max(which(up.sub.mat.col.sums==min(up.sub.mat.col.sums)))+(dim(mod.mat)[1]-length(up.sub.mat.col.sums))
          lower.bound<-min(which(up.sub.mat.col.sums==min(up.sub.mat.col.sums)))+(dim(mod.mat)[1]-length(up.sub.mat.col.sums))
          if(upper.bound==dim(mat)[1])
          {
            output[6]<-"upper.erad"
            output[7]<-virulence.steps[lower.bound]
          }
        }
        
        low.sub.mat.cols<-which(virulence.steps<output[2])
        low.sub.mat<-mod.mat[low.sub.mat.cols,low.sub.mat.cols]
        low.sub.mat.col.sums<-colSums(low.sub.mat)
        
        if (length(which(low.sub.mat.col.sums==min(low.sub.mat.col.sums)))>=2)
        {
          upper.bound<-max(which(low.sub.mat.col.sums==min(low.sub.mat.col.sums)))
          lower.bound<-min(which(low.sub.mat.col.sums==min(low.sub.mat.col.sums)))
          
          if(lower.bound==1)
          {
            output[8]<-"lower.erad"
            output[9]<-virulence.steps[upper.bound]
          }
        }
      }
    }
  }
  
  if (is.na(output[5]) && is.na(col.max) && is.na(col.min) && is.na(row.max) && is.na(row.min)) {output[5]<-"global eradication"}
  
  return(output)
}





