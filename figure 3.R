# setup
source("~/Documents/GitHub/covid_vaccines_virulence_evolution/functions.R")

## plotting function for visualizing results matrix
sub.plot.func<-function(matrix,cols,col.vals)
{
  plot(0,0,type="n",xlim=c(-(1/(res-1))/2,1+(1/(res-1))/2),ylim=c(-(1/(res-1))/2,1+(1/(res-1))/2),xlab=expression('r'[U]),ylab=expression('r'[L]),cex.lab=2,axes=F)
  axis(1,at=seq(0,1,.2),labels = seq(.5,1,.1))
  axis(2,at=seq(0,1,.2),labels = seq(.5,1,.1))
  box()
  xx<-seq(0,1,length.out = res)
  yy<-seq(0,1,length.out = res)
  for(i in 1:res)
  {
    for(j in 1:res)
    {
      rect(xx[i]-(1/(res-1))/2,yy[j]-(1/(res-1))/2,xx[i]+(1/(res-1))/2,yy[j]+(1/(res-1))/2,col = cols[which.min(abs(matrix[i,j]-col.vals))],border=NA)
    }
  }
}

plot.func<-function(plot.data.,colors.=colors,col.vals.=col.vals,optim.vir)
{
  alpha.ess.mat.data<-as.numeric(unlist(plot.data.$alpha.ess))
  alpha.ess.mat.data[which(plot.data.$pip.motif=="selection for hypervirulence")]<-(-2)
  alpha.ess.mat.data[which(plot.data.$pip.motif=="global eradication")]<-(-1)
  alpha.ess.mat<-matrix(alpha.ess.mat.data,res,res,byrow = T)
  sub.plot.func(alpha.ess.mat,colors.,col.vals.)
  
  epi.result.mat<-matrix(NA,res,res)
  sub.mat<-data.frame(x=c(),y=c())
  for (i in 1:nrow(plot.data))
  {
    xindex<-which(rUv.steps==plot.data[i,'rUv'])
    yindex<-which(rLv.steps==plot.data[i,'rLv'])
    xx<-seq(0,1,length.out = res)
    yy<-seq(0,1,length.out = res)
    if(plot.data[i,"Re.alpha.delta.start"] >= 1 & (isTRUE(plot.data[i,"alpha.ess"] > optim.vir) | (plot.data[i,"pip.motif"]=="selection for hypervirulence"))) {outcome<-1}
    if(plot.data[i,"Re.alpha.delta.start"] >= 1 & isTRUE(plot.data[i,"alpha.ess"] <= optim.vir)) {outcome<-2}
    if(plot.data[i,"Re.alpha.delta.start"] < 1 & (isTRUE(plot.data[i,"alpha.ess"] > optim.vir) | (plot.data[i,"pip.motif"]=="selection for hypervirulence"))) {outcome<-3}
    if(plot.data[i,"Re.alpha.delta.start"] < 1 & isTRUE(plot.data[i,"alpha.ess"] <= optim.vir)) {outcome<-4}
    if(plot.data[i,"Re.alpha.delta.start"] < 1 & plot.data[i,"pip.motif"]=="global eradication") {outcome<-5}
    epi.result.mat[xindex,yindex]<-outcome
    if(outcome==3) {sub.mat<-rbind(sub.mat,data.frame(x=xx[xindex],y=yy[yindex]))}
  }
  sub.mat<-round(sub.mat,4)
  poly.cords<-data.frame(x=c(),y=c())
  for(i in unique(sub.mat$x))
  {
    sub.col<-sub.mat[which(sub.mat$x==i),]
    poly.cords<-rbind(poly.cords,data.frame(x=i-(1/(res-1))/2,y=min(sub.col$y)-(1/(res-1))/2))
    poly.cords<-rbind(poly.cords,data.frame(x=i-(1/(res-1))/2,y=max(sub.col$y)+(1/(res-1))/2))
    poly.cords<-rbind(poly.cords,data.frame(x=i+(1/(res-1))/2,y=min(sub.col$y)-(1/(res-1))/2))
    poly.cords<-rbind(poly.cords,data.frame(x=i+(1/(res-1))/2,y=max(sub.col$y)+(1/(res-1))/2))
  }
  
  if(nrow(poly.cords)>0)
  {
    poly.cords<-poly.cords[order(poly.cords$y,poly.cords$x),]
    upper<-poly.cords[1:(nrow(poly.cords)/2),]
    lower<-poly.cords[(1+nrow(poly.cords)/2):nrow(poly.cords),]
    lower<-lower[nrow(lower):1,]
    poly.cords.sorted<-rbind(upper,lower)
    
    polygon(poly.cords.sorted$x,poly.cords.sorted$y,lwd=1,density = 10)
  }
}



## params
A.plot<-c(-2,-1,seq(0,.2,.00025))
colors<-c("grey","white",rev(magma(length(A.plot))))
col.vals<-A.plot

res<-11
rUv.steps<-rLv.steps<-seq(.5,1,length.out = res)

# plot
#layout(matrix(c(1,1,2,2,3,3,10,4,4,5,5,6,6,10,7,7,8,8,9,9,11),3,7,byrow = T))
layout(matrix(c(1,1,2,2,3,3,10,1,1,2,2,3,3,10,4,4,5,5,6,6,10,4,4,5,5,6,6,11,7,7,8,8,9,9,11,7,7,8,8,9,9,11),6,7,byrow=T))
par(mar=c(2,2,2,2),oma=c(4,8,4,0))

## alpha optim = 0.00875 (intermediate between alpha and delta strain virulence)

## 50% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.5alpha.optim0.00875.RDS")
plot.func(plot.data,optim.vir = .00875)

mtext("50% vaccinated",line=2,cex=1.25)
mtext(expression(alpha['optim']*' = '*(alpha[italic("delta")]+alpha[italic("alpha")])/2),side=2,line=7,cex=1.25)

## 75% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.75alpha.optim0.00875.RDS")
plot.func(plot.data,optim.vir = .00875)

mtext("75% vaccinated",line=2,cex=1.25)

## 90% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.9alpha.optim0.00875.RDS")
plot.func(plot.data,optim.vir = .00875)

mtext("90% vaccinated",line=2,cex=1.25)

## alpha optim = 0.01 (delta strain virulence)

## 50% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.5alpha.optim0.01.RDS")
plot.func(plot.data,optim.vir = .01)

mtext(expression(alpha['optim']*' = '*alpha[italic("delta")]),side=2,line=7,cex=1.25)

## 75% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.75alpha.optim0.01.RDS")
plot.func(plot.data,optim.vir = .01)
## 90% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.9alpha.optim0.01.RDS")
plot.func(plot.data,optim.vir = .01)

## alpha optim = 0.02 (twice as virulent as delta)

## 50% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.5alpha.optim0.02.RDS")
plot.func(plot.data,optim.vir = .02)
mtext(expression(alpha['optim']*' = 2*'*alpha[italic("delta")]),side=2,line=7,cex=1.25)
## 75% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.75alpha.optim0.02.RDS")
plot.func(plot.data,optim.vir = .02)
## 90% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.9alpha.optim0.02.RDS")
plot.func(plot.data,optim.vir = .02)


mtext(expression('lower respiratory tract protection (r'["L,V"]*')'),side = 2,line=1.5,cex=1.5,outer=T)
mtext(expression('upper respiratory tract protection (r'["U,V"]*')'),side = 1,line=2,cex=1.5,outer=T,adj=2.8/7)

# legend

par(mar=c(0,1,2,5))
yy<-seq(0,length(col.vals)-2,1)
plot(0,0,type="n",xlim=c(0,1),ylim=c(.5,length(col.vals)-2+.5),xlab="",ylab="",axes=F)
color.legend(0,0,1,length(col.vals)-2,legend=NULL,colors,gradient="y")
axis(4,at=seq(0,length(col.vals)-2,length.out = 5),labels = col.vals[seq(3,length(col.vals),length.out = 5)])
mtext(expression(alpha["ESS"]),side=3,line =0,cex=2)

par(mar=c(0,1,5,5),xpd=T)
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=F)
polygon(c(0,0,1,1),c(0,.2,.2,0),density=10)
text(.5,.2,expression(atop(R[E[italic("delta")]]<1,"\nbefore evolution")),pos=3,cex=1.25)
rect(0,.4,1,.6,col="white")
text(.5,.6,"evolution proof\nherd immunity",pos=3,cex=1.25)
rect(0,.8,1,1,col="grey")
text(.5,1,"unbounded\nselection for\nincreased virulence",pos=3,cex=1.25)
par(xpd=F)
## copy to clipboard with width of 1285, height of 868

