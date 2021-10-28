# setup
col.trans<-.6

colors1<-rev(magma(length(A)))
col.vals1<-A

res<-11
rUv.steps<-rLv.steps<-seq(.5,1,length.out = res)

# plot
par(mfrow=c(3,4))

## alpha optim = 0.00875 (intermediate between alpha and delta strain virulence)

## 50% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega10p.vacc0.5alpha.optim0.00875.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)
## 75% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega10p.vacc0.75alpha.optim0.00875.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)
## 90% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega10p.vacc0.9alpha.optim0.00875.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)
## 99% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega10p.vacc0.99alpha.optim0.00875.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)

## alpha optim = 0.01 (delta strain virulence)

## 50% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega10p.vacc0.5alpha.optim0.01.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)
## 75% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega10p.vacc0.75alpha.optim0.01.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)
## 90% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega10p.vacc0.9alpha.optim0.01.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)
## 99% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega10p.vacc0.99alpha.optim0.01.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)

## alpha optim = 0.02 (twice as virulent as delta)

## 50% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega10p.vacc0.5alpha.optim0.02.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)
## 75% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega10p.vacc0.75alpha.optim0.02.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)
## 90% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega10p.vacc0.9alpha.optim0.02.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)
## 99% vaccinated
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega10p.vacc0.99alpha.optim0.02.RDS")
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)





mtext(expression('lower respiratory tract protection (r'["L,V"]*')'),side = 2,line=1.5,cex=1.5,outer=T)
mtext(expression('upper respiratory tract protection (r'["U,V"]*')'),side = 1,line=2,cex=1.5,outer=T,adj=3/7)

## copy to clipboard with width of 1051, height of 843


