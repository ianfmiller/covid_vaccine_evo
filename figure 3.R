# setup
col.trans<-.6

colors1<-rev(magma(length(virulence.steps)))
col.vals1<-virulence.steps

colors2<-c("yellow","blue","green","darkgreen")
col.vals2<-c(1,2,3,4)

s.blues3<-hsv(.666,1,1,seq(1,.001,length.out = 200)^col.trans)
s.reds3<-hsv(1,1,1,seq(.001,1,length.out = 200)^col.trans)
s.colors3<-c(s.blues3,"white",s.reds3)
s.col.vals3<-seq(-7.5,7.5,length.out = 401)

s.colors<-list(s.colors1,s.colors2,s.colors3)
s.col.vals<-list(s.col.vals1,s.col.vals2,s.col.vals3)

res<-21
rUv.steps<-rLv.steps<-seq(.5,1,length.out = res)

# plot

## alpha optim = 0.00875 (intermediate between alpha and delta strain virulence)

par(mfrow=c(2,2))
plot.data<-readRDS("~/Documents/GitHub/covid_vaccines_virulence_evolution/sim.data/omega10p.vacc0.99alpha.optim0.02.RDS")
unique(plot.data$pip.motif)
hist(plot.data$Re.alpha.delta.start,xlim=c(0,5))
alpha.ess.mat<-matrix(as.numeric(unlist(plot.data$alpha.ess)),res,res,byrow = T)
plot.result(alpha.ess.mat,colors1,col.vals1)

epi.result.mat<-matrix(NA,res,res)
for (i in 1:nrow(plot.data))
{
  xindex<-which(rUv.steps==plot.data[i,'rUv'])
  yindex<-which(rLv.steps==plot.data[i,'rLv'])
  if(plot.data[i,"Re.alpha.delta.start"] >= 1 & plot.data[i,"alpha.ess"] > .00875) {outcome<-1}
  if(plot.data[i,"Re.alpha.delta.start"] >= 1 & plot.data[i,"alpha.ess"] <= .00875) {outcome<-2}
  if(plot.data[i,"Re.alpha.delta.start"] < 1 & plot.data[i,"alpha.ess"] > .00875) {outcome<-3}
  if(plot.data[i,"Re.alpha.delta.start"] < 1 & plot.data[i,"alpha.ess"] <= .00875) {outcome<-4}
  
  epi.result.mat[xindex,yindex]<-outcome
}

plot.result(epi.result.mat,colors2,col.vals2)

## 50% vaccinated

## 75% vaccinated

## 90% vaccinated

## alpha optim = 0.01 (delta strain virulence)

## 50% vaccinated

## 75% vaccinated

## 90% vaccinated

## alpha optim = 0.02 (twice as virulent as delta)

## 50% vaccinated

## 75% vaccinated

## 90% vaccinated

plot.mat.R0.obs<-matrix(data2$Re.alpha.delta.epi.equi,res,res,byrow = T) #populate matricies
plot.s(plot.mat.R0.obs,s.colors[[color.index]],s.col.vals[[color.index]])
contour(plot.mat.R0.mutant-plot.mat.R0.obs,add=T)
mtext(expression('r'[L]),side = 2,line=2.5)
mtext("10% vaccinated",line=2,cex=1.25)

plot.mat.alpha.ess<-matrix(,res,res,byrow = T) #populate matricies
plot.mat.R0.mutant<-matrix(R0.mutant.vec,res,res,byrow = T) #populate matricies
s.mat<-plot.mat.R0.mutant-plot.mat.R0.obs
plot.s(s.mat,s.colors[[color.index]],s.col.vals[[color.index]])
contour(plot.mat.R0.mutant-plot.mat.R0.obs,add=T)
mtext("50% vaccinated",line=2,cex=1.25)
#mtext(expression('selection for '*alpha*' = 0.01'),side=3,line=4,font=2,cex=1.2)



  
plot.mat.R0.obs<-matrix(R0.obs.vec,res,res,byrow = T) #populate matricies
plot.mat.R0.mutant<-matrix(R0.mutant.vec,res,res,byrow = T) #populate matricies
s.mat<-plot.mat.R0.mutant-plot.mat.R0.obs
plot.s(s.mat,s.colors[[color.index]],s.col.vals[[color.index]])
contour(plot.mat.R0.mutant-plot.mat.R0.obs,add=T)
mtext("90% vaccinated",line=2,cex=1.25)





plot.mat.R0.obs<-matrix(R0.obs.vec,res,res,byrow = T) #populate matricies
plot.mat.R0.mutant<-matrix(R0.mutant.vec,res,res,byrow = T) #populate matricies
s.mat<-plot.mat.R0.mutant-plot.mat.R0.obs
plot.s(s.mat,s.colors[[color.index]],s.col.vals[[color.index]])
contour(plot.mat.R0.mutant-plot.mat.R0.obs,add=T)
#mtext(expression('r'[L]),side = 2,line=2.5)
mtext(expression(alpha['optim']*' = '*alpha[B.1.1.7]),side=2,line=7,cex=1.25)



plot.mat.R0.obs<-matrix(R0.obs.vec,res,res,byrow = T) #populate matricies
plot.mat.R0.mutant<-matrix(R0.mutant.vec,res,res,byrow = T) #populate matricies
s.mat<-plot.mat.R0.mutant-plot.mat.R0.obs
plot.s(s.mat,s.colors[[color.index]],s.col.vals[[color.index]])
contour(plot.mat.R0.mutant-plot.mat.R0.obs,add=T)
#mtext(expression('selection for '*alpha*' = 0.01'),side=3,line=4,font=2,cex=1.2)



plot.mat.R0.obs<-matrix(R0.obs.vec,res,res,byrow = T) #populate matricies
plot.mat.R0.mutant<-matrix(R0.mutant.vec,res,res,byrow = T) #populate matricies
s.mat<-plot.mat.R0.mutant-plot.mat.R0.obs
plot.s(s.mat,s.colors[[color.index]],s.col.vals[[color.index]])
contour(plot.mat.R0.mutant-plot.mat.R0.obs,add=T)



plot.mat.R0.obs<-matrix(R0.obs.vec,res,res,byrow = T) #populate matricies
plot.mat.R0.mutant<-matrix(R0.mutant.vec,res,res,byrow = T) #populate matricies
s.mat<-plot.mat.R0.mutant-plot.mat.R0.obs
plot.s(s.mat,s.colors[[color.index]],s.col.vals[[color.index]])
contour(plot.mat.R0.mutant-plot.mat.R0.obs,add=T)
mtext(expression(alpha['optim']*' = 1.25*'*alpha[ansc]),side=2,line=7,cex=1.25)



plot.mat.R0.obs<-matrix(R0.obs.vec,res,res,byrow = T) #populate matricies
plot.mat.R0.mutant<-matrix(R0.mutant.vec,res,res,byrow = T) #populate matricies
s.mat<-plot.mat.R0.mutant-plot.mat.R0.obs
plot.s(s.mat,s.colors[[color.index]],s.col.vals[[color.index]])
contour(plot.mat.R0.mutant-plot.mat.R0.obs,add=T)



plot.mat.R0.obs<-matrix(R0.obs.vec,res,res,byrow = T) #populate matricies
plot.mat.R0.mutant<-matrix(R0.mutant.vec,res,res,byrow = T) #populate matricies
s.mat<-plot.mat.R0.mutant-plot.mat.R0.obs
plot.s(s.mat,s.colors[[color.index]],s.col.vals[[color.index]])
contour(plot.mat.R0.mutant-plot.mat.R0.obs,add=T)




mtext(expression('lower respiratory tract protection (r'["L,V"]*')'),side = 2,line=1.5,cex=1.5,outer=T)
mtext(expression('upper respiratory tract protection (r'["U,V"]*')'),side = 1,line=2,cex=1.5,outer=T,adj=3/7)

## copy to clipboard with width of 1051, height of 843


