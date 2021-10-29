# setup
source("~/Documents/GitHub/covid_vaccines_virulence_evolution/functions.R")

A.plot<-c(-2,-1,seq(0,.2,.00025))
colors<-c("grey","white",rev(magma(length(A.plot))))
col.vals<-A.plot

res<-11

# plot
layout(matrix(c(1,1,2,2,3,3,10,4,4,5,5,6,6,10,7,7,8,8,9,9,11),3,7,byrow = T))
par(mar=c(2,2,2,2),oma=c(4,8,4,0))

## alpha optim = 0.00875 (intermediate between alpha and delta strain virulence)

## 50% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.5alpha.optim0.00875.RDS")
alpha.ess.mat.data<-as.numeric(unlist(plot.data$alpha.ess))
alpha.ess.mat.data[which(plot.data$pip.motif=="selection for hypervirulence")]<-(-2)
alpha.ess.mat.data[which(plot.data$pip.motif=="global eradication")]<-(-1)
alpha.ess.mat<-matrix(alpha.ess.mat.data,res,res,byrow = T)
plot.result(alpha.ess.mat,colors,col.vals)

mtext("50% vaccinated",line=2,cex=1.25)
mtext(expression(alpha['optim']*' = '*(alpha["deltaa"]+alpha["alpha"])/2),side=2,line=7,cex=1.25)

## 75% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.75alpha.optim0.00875.RDS")
alpha.ess.mat.data<-as.numeric(unlist(plot.data$alpha.ess))
alpha.ess.mat.data[which(plot.data$pip.motif=="selection for hypervirulence")]<-(-2)
alpha.ess.mat.data[which(plot.data$pip.motif=="global eradication")]<-(-1)
alpha.ess.mat<-matrix(alpha.ess.mat.data,res,res,byrow = T)
plot.result(alpha.ess.mat,colors,col.vals)

mtext("75% vaccinated",line=2,cex=1.25)

## 90% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.9alpha.optim0.00875.RDS")
alpha.ess.mat.data<-as.numeric(unlist(plot.data$alpha.ess))
alpha.ess.mat.data[which(plot.data$pip.motif=="selection for hypervirulence")]<-(-2)
alpha.ess.mat.data[which(plot.data$pip.motif=="global eradication")]<-(-1)
alpha.ess.mat<-matrix(alpha.ess.mat.data,res,res,byrow = T)
plot.result(alpha.ess.mat,colors,col.vals)

mtext("90% vaccinated",line=2,cex=1.25)

## alpha optim = 0.01 (delta strain virulence)

## 50% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.5alpha.optim0.01.RDS")
alpha.ess.mat.data<-as.numeric(unlist(plot.data$alpha.ess))
alpha.ess.mat.data[which(plot.data$pip.motif=="selection for hypervirulence")]<-(-2)
alpha.ess.mat.data[which(plot.data$pip.motif=="global eradication")]<-(-1)
alpha.ess.mat<-matrix(alpha.ess.mat.data,res,res,byrow = T)
plot.result(alpha.ess.mat,colors,col.vals)

mtext(expression(alpha['optim']*' = '*alpha["deltaa"]),side=2,line=7,cex=1.25)

## 75% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.75alpha.optim0.01.RDS")
alpha.ess.mat.data<-as.numeric(unlist(plot.data$alpha.ess))
alpha.ess.mat.data[which(plot.data$pip.motif=="selection for hypervirulence")]<-(-2)
alpha.ess.mat.data[which(plot.data$pip.motif=="global eradication")]<-(-1)
alpha.ess.mat<-matrix(alpha.ess.mat.data,res,res,byrow = T)
plot.result(alpha.ess.mat,colors,col.vals)
## 90% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.9alpha.optim0.01.RDS")
alpha.ess.mat.data<-as.numeric(unlist(plot.data$alpha.ess))
alpha.ess.mat.data[which(plot.data$pip.motif=="selection for hypervirulence")]<-(-2)
alpha.ess.mat.data[which(plot.data$pip.motif=="global eradication")]<-(-1)
alpha.ess.mat<-matrix(alpha.ess.mat.data,res,res,byrow = T)
plot.result(alpha.ess.mat,colors,col.vals)

## alpha optim = 0.02 (twice as virulent as delta)

## 50% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.5alpha.optim0.02.RDS")
alpha.ess.mat.data<-as.numeric(unlist(plot.data$alpha.ess))
alpha.ess.mat.data[which(plot.data$pip.motif=="selection for hypervirulence")]<-(-2)
alpha.ess.mat.data[which(plot.data$pip.motif=="global eradication")]<-(-1)
alpha.ess.mat<-matrix(alpha.ess.mat.data,res,res,byrow = T)
plot.result(alpha.ess.mat,colors,col.vals)
mtext(expression(alpha['optim']*' = 2*'*alpha["deltaa"]),side=2,line=7,cex=1.25)
## 75% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.75alpha.optim0.02.RDS")
alpha.ess.mat.data<-as.numeric(unlist(plot.data$alpha.ess))
alpha.ess.mat.data[which(plot.data$pip.motif=="selection for hypervirulence")]<-(-2)
alpha.ess.mat.data[which(plot.data$pip.motif=="global eradication")]<-(-1)
alpha.ess.mat<-matrix(alpha.ess.mat.data,res,res,byrow = T)
plot.result(alpha.ess.mat,colors,col.vals)
## 90% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.9alpha.optim0.02.RDS")
alpha.ess.mat.data<-as.numeric(unlist(plot.data$alpha.ess))
alpha.ess.mat.data[which(plot.data$pip.motif=="selection for hypervirulence")]<-(-2)
alpha.ess.mat.data[which(plot.data$pip.motif=="global eradication")]<-(-1)
alpha.ess.mat<-matrix(alpha.ess.mat.data,res,res,byrow = T)
plot.result(alpha.ess.mat,colors,col.vals)


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
rect(0,.2,1,.4,col="white")
text(.5,.4,"evolution proof\nherd immunity",pos=3,cex=1.25)
rect(0,.7,1,.9,col="grey")
text(.5,.9,"unbounded\nselection for\nincreased virulence",pos=3,cex=1.25)
par(xpd=F)
## copy to clipboard with width of 1051, height of 843

