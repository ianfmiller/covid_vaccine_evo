# setup
A.plot<-c(-2,-1,seq(.015,.2,.00025))
colors1<-c("black","grey",rev(magma(length(A.plot))))
col.vals1<-A.plot

res<-11
rUv.steps<-rLv.steps<-seq(.5,1,length.out = res)

# plot
par(mfrow=c(3,4),oma=c(6,6,0,0))

## alpha optim = 0.00875 (intermediate between alpha and delta strain virulence)

## 50% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.5alpha.optim0.00875.RDS")
alpha.ess.mat.data<-as.numeric(unlist(plot.data$alpha.ess))
alpha.ess.mat.data[which(plot.data$pip.motif=="selection for hypervirulence")]<-(-2)
alpha.ess.mat.data[which(plot.data$pip.motif=="global eradication")]<-(-1)
plot.result(alpha.ess.mat,colors1,col.vals1)
## 75% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.75alpha.optim0.00875.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)
## 90% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.9alpha.optim0.00875.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)
## 99% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.99alpha.optim0.00875.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)

## alpha optim = 0.01 (delta strain virulence)

## 50% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.5alpha.optim0.01.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)
## 75% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.75alpha.optim0.01.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)
## 90% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.9alpha.optim0.01.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)
## 99% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.99alpha.optim0.01.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)

## alpha optim = 0.02 (twice as virulent as delta)

## 50% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.5alpha.optim0.02.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)
## 75% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.75alpha.optim0.02.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)
## 90% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.9alpha.optim0.02.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)
## 99% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega0p.vacc0.99alpha.optim0.02.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)





mtext(expression('lower respiratory tract protection (r'["L,V"]*')'),side = 2,line=1.5,cex=1.5,outer=T)
mtext(expression('upper respiratory tract protection (r'["U,V"]*')'),side = 1,line=2,cex=1.5,outer=T,adj=4/8)

## copy to clipboard with width of 1051, height of 843


