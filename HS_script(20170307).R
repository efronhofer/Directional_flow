#..........................................#
#...........HetScape project...............#
#..........................................#

#..........................................................................................................................................#
#... Collaborators: Eric Harvey, Isabelle Gounand, Chelsea Little, Emanuel Fronhofer, Florian Altermatt                  #
#... Author of the script: Eric Harvey - contact: eric.harvey@eawag.ch                                                                                                #
#... Date creation: 2016-09-14                                                                                                           #                                                                       #
#..........................................................................................................................................#

#... Clear variables from environment
rm(list=ls())
search()
#detach(pos=2)

#... Set directories paths
datapath = "~/Documents/Research/Eawag/Projects/5.HetScape/3.Analysis/"
graphpath = "~/Documents/Research/Eawag/Projects/5.HetScape/4.Results/"

#...Load needed packages
library(sciplot)
library(car)
library(MASS)
library(plyr)
library(coefplot2)
library(nlme)
library(scales)
library(vegan)
#citation(package = "base", lib.loc = NULL)

#... Load the data
HS.data = read.delim(paste0(datapath,"HS_alldata_20161031.txt"))

#...Create new variable for experimental days
HS.data$Day = factor(c(rep("0",160),rep("3",160),rep("7",160),rep("10",160),rep("14",160),rep("17",160),rep("21",160)), 
                     labels = c("0","3","7","10","14","17","21"))

#########################################################################################################
#...Data structure

#......An uncontrolled event happened after day 14 that left a strong signature on bacteria density and oxygen concentration, therefore we decided to exclude these two last days from further analyses
with(HS.data,lineplot.CI(Day,log(Bact_per_ul),xlab="Days",ylab="Bacteria density"))
with(HS.data,lineplot.CI(Day,Oxygen,xlab="Days",ylab="Oxygen (%)"))

#.....Remove Day 0 and the two last days for further analyses
sel.day = c("3","7","10","14")
HS.data2 = HS.data[HS.data$Day %in% sel.day,]
HS.data2 = droplevels(HS.data2)
#.....Convert Day into a numerical variable
HS.data2$Day = as.numeric(levels(HS.data2$Day))[HS.data2$Day] 

#.....Choose dilution to keep for further analyses (0 or 0.33) 
HS.data2$Dilution = as.factor(HS.data2$Dilution)
sel.dil = c("0")
HS.data3 = HS.data2[HS.data2$Dilution %in% sel.dil,]
HS.data3 = droplevels(HS.data3)

#.....Seperate Ecosystem 1 and 2 (Upstream and downstream respectively)
Eco1.data = subset(HS.data3,Ecosystem == 1)
Eco1.data = droplevels(Eco1.data)
#rename factor levels to correspond to manuscript
Eco1.data$Interaction = with(Eco1.data,revalue(Interaction, c("C.alone"="Monoculture", "B.alone"="Resource.alone")))
Eco1.data$Interaction = factor(Eco1.data$Interaction,levels=c("Monoculture","Competition","Facilitation","Predation","Resource.alone"))
Eco2.data = subset(HS.data3,Ecosystem == 2)
Eco2.data = droplevels(Eco2.data)
Eco2.data$Interaction = with(Eco2.data,revalue(Interaction, c("C.alone"="Monoculture", "B.alone"="Resource.alone")))
Eco2.data$Interaction = factor(Eco2.data$Interaction,levels=c("Monoculture","Competition","Facilitation","Predation","Resource.alone"))

#....Separate Ecosystem 2 with the consumer (Tet) present 
Eco2.data.tet = subset(Eco2.data,Composition == "Tet")
Eco2.data.tet = droplevels(Eco2.data.tet)

#########################################################################################################
#...Statistical Analysis for Ecosystem 2 (focus of this experiment)

#..Identify probability distribution for Bacteria density in Eco 2
summary(Eco2.data$Bact_per_ul)
with(Eco2.data,plot(density(Bact_per_ul)))
with(Eco2.data,plot(density(decostand(Bact_per_ul,"log"))))
with(Eco2.data,plot(density(sqrt(Bact_per_ul))))
with(Eco2.data,qqp(Bact_per_ul,"norm"))
with(Eco2.data,qqp(Bact_per_ul,"lnorm"))
poisson <- with(Eco2.data,fitdistr(Bact_per_ul, "Poisson"))
with(Eco2.data,qqp(Bact_per_ul, "pois",poisson$estimate))

#..Identify probability distribution for Tetrahymena density in Eco 2
summary(Eco2.data.tet$Tet)
with(Eco2.data.tet,plot(density(Tet)))
with(Eco2.data.tet,plot(density(log(Tet))))
with(Eco2.data.tet,qqp(Tet,"norm"))
with(Eco2.data.tet,qqp(Tet,"lnorm"))
poisson <- with(Eco2.data.tet,fitdistr(Tet, "Poisson"))
with(Eco2.data.tet,qqp(Tet, "pois",poisson$estimate))
gamma = with(Eco2.data.tet,fitdistr(Tet, "gamma")) #only for value larger than zero
with(Eco2.data.tet,qqp(Tet, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]]))

#..Statistical models for Bacteria density in Ecosystem 2

#Linear mixed model (AIC-based simplification procedure - see Methods)

mod1.3 = with(Eco2.data,lme(log(Bact_per_ul) ~ Composition*Interaction*Day,random = ~ Day|Replicate,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim")))
summary(mod1.3)$tTable
summary(mod1.3)$AIC

mod1.3.1 = with(Eco2.data,update(mod1.3, ~. - Composition:Interaction:Day))
summary(mod1.3.1)$AIC
anova(mod1.3,mod1.3.1)

mod1.3.2 = with(Eco2.data,update(mod1.3.1, ~. - Interaction:Day))
anova(mod1.3.1,mod1.3.2)
summary(mod1.3.2)$AIC

mod1.3.3 = with(Eco2.data,update(mod1.3.2, ~. - Composition:Day))
anova(mod1.3.2,mod1.3.3)
summary(mod1.3.3)$AIC

mod1.3.4 = with(Eco2.data,update(mod1.3.3, ~. - Composition:Interaction))
anova(mod1.3.3,mod1.3.4)
summary(mod1.3.4)$AIC

mod1.3.5 = with(Eco2.data,update(mod1.3.3, ~. - Day))
anova(mod1.3.3,mod1.3.5)
summary(mod1.3.5)$AIC
summary(mod1.3.5)$tTable
coefplot2(mod1.3.5,vertical=T)

#Rerun final model with REML estimator
mod1.4 = with(Eco2.data,lme(log(Bact_per_ul) ~ Composition*Interaction + Day,random = ~ Day|Replicate,method="REML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim")))
summary(mod1.4)$tTable
coefplot2(mod1.4,vertical=T)
dev.off()
intervals(mod1.4)

#Extract and plot model prediction + CI (http://stackoverflow.com/questions/14358811/extract-prediction-band-from-lme-fit)
X= rep(seq(3,14,len=4),times=1,each=40)
new.dat = data.frame(Day=X,Composition= Eco2.data$Composition,Interaction=Eco2.data$Interaction)
new.dat$pred = predict(mod1.4,newdata=new.dat,level=0)
Designmat = model.matrix(eval(eval(mod1.4$call$fixed)[-2]),new.dat[-ncol(new.dat)])
predvar = diag(Designmat %*% mod1.4$varFix %*% t(Designmat))
new.dat$SE = sqrt(predvar)

#Figure 2

colo.s = c("#99999930", "#E69F0030", "#56B4E930", "#009E7330","#CC79A730")
colo.l = c("#999999", "#E69F00", "#56B4E9", "#009E73","#CC79A7")
#Figure frame
with(Eco2.data,plot(Day,log(Bact_per_ul),type="n",xaxt="n",yaxt="n"))
axis(side = 2, at =c(10.13,11.44,12,12.5,13.12),labels=c(round(exp(10.13),-4),round(exp(11.44),-4),round(exp(12),-4),round(exp(12.5),-4),round(exp(13.12),-3)))
axis(side = 1, at =c(0,3,7,10,14))

#Without Tet

#Draw Shadings first and then lines (to avoid having shadings hiding lines)
x1 = new.dat$Day[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Monoculture")]
y1.1 = new.dat$pred[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Monoculture")]+new.dat$SE[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Monoculture")]
y2.1 = new.dat$pred[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Monoculture")]-new.dat$SE[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Monoculture")]
polygon(c(x1,rev(x1)),c(y2.1,rev(y1.1)),col=colo.s[2],border=0)

x2 = new.dat$Day[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Resource.alone")]
y1.2 = new.dat$pred[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Resource.alone")]+new.dat$SE[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Resource.alone")]
y2.2 = new.dat$pred[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Resource.alone")]-new.dat$SE[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Resource.alone")]
polygon(c(x2,rev(x2)),c(y2.2,rev(y1.2)),col=colo.s[1],border=0)

x3 =  new.dat$Day[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Competition")]
y1.3 = new.dat$pred[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Competition")]+new.dat$SE[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Competition")]
y2.3 = new.dat$pred[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Competition")]-new.dat$SE[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Competition")]
polygon(c(x3,rev(x3)),c(y2.3,rev(y1.3)),col=colo.s[3],border=0)

x4 = new.dat$Day[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Facilitation")]
y1.4 = new.dat$pred[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Facilitation")]+new.dat$SE[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Facilitation")]
y2.4 = new.dat$pred[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Facilitation")]-new.dat$SE[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Facilitation")]
polygon(c(x4,rev(x4)),c(y2.4,rev(y1.4)),col=colo.s[4],border=0)

x5 = new.dat$Day[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Predation")]
y1.5 = new.dat$pred[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Predation")]+new.dat$SE[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Predation")]
y2.5 = new.dat$pred[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Predation")]-new.dat$SE[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Predation")]
polygon(c(x5,rev(x5)),c(y2.5,rev(y1.5)),col=colo.s[5],border=0)

lines(new.dat$Day[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Monoculture")],new.dat$pred[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Monoculture")],lwd=2,col=colo.l[2])
points(new.dat$Day[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Monoculture")],log(Eco2.data$Bact_per_ul[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Monoculture")]),lwd=2,col=alpha(colo.l[2],0.4),pch=16)
lines(new.dat$Day[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Resource.alone")],new.dat$pred[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Resource.alone")],lwd=2,col=colo.l[1])
points(new.dat$Day[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Resource.alone")],log(Eco2.data$Bact_per_ul[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Resource.alone")]),lwd=2,col=alpha(colo.l[1],0.4),pch=16)
lines(new.dat$Day[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Competition")],new.dat$pred[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Competition")],lwd=2,col=colo.l[3])
points(new.dat$Day[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Competition")],log(Eco2.data$Bact_per_ul[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Competition")]),lwd=2,col=alpha(colo.l[3],0.4),pch=16)
lines(new.dat$Day[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Facilitation")],new.dat$pred[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Facilitation")],lwd=2,col=colo.l[4])
points(new.dat$Day[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Facilitation")],log(Eco2.data$Bact_per_ul[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Facilitation")]),lwd=2,col=alpha(colo.l[4],0.4),pch=16)
lines(new.dat$Day[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Predation")],new.dat$pred[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Predation")],lwd=2,col=colo.l[5])
points(new.dat$Day[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Predation")],log(Eco2.data$Bact_per_ul[which(new.dat$Composition == "Bact" & new.dat$Interaction=="Predation")]),lwd=2,col=alpha(colo.l[5],0.4),pch=16)


#With Tet
x6 = new.dat$Day[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Monoculture")]
y1.6 = new.dat$pred[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Monoculture")]+new.dat$SE[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Monoculture")]
y2.6 = new.dat$pred[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Monoculture")]-new.dat$SE[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Monoculture")]
polygon(c(x6,rev(x6)),c(y2.6,rev(y1.6)),col=colo.s[2],border=0)

x7 = new.dat$Day[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Resource.alone")]
y1.7 = new.dat$pred[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Resource.alone")]+new.dat$SE[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Resource.alone")]
y2.7 = new.dat$pred[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Resource.alone")]-new.dat$SE[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Resource.alone")]
polygon(c(x7,rev(x7)),c(y2.7,rev(y1.7)),col=colo.s[1],border=0)

x8 = new.dat$Day[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Competition")]
y1.8 = new.dat$pred[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Competition")]+new.dat$SE[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Competition")]
y2.8 = new.dat$pred[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Competition")]-new.dat$SE[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Competition")]
polygon(c(x8,rev(x8)),c(y2.8,rev(y1.8)),col=colo.s[3],border=0)

x9 = new.dat$Day[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Facilitation")]
y1.9 = new.dat$pred[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Facilitation")]+new.dat$SE[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Facilitation")]
y2.9 = new.dat$pred[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Facilitation")]-new.dat$SE[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Facilitation")]
polygon(c(x9,rev(x9)),c(y2.9,rev(y1.9)),col=colo.s[4],border=0)

x10 = new.dat$Day[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Predation")]
y1.10 = new.dat$pred[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Predation")]+new.dat$SE[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Predation")]
y2.10 = new.dat$pred[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Predation")]-new.dat$SE[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Predation")]
polygon(c(x10,rev(x10)),c(y2.10,rev(y1.10)),col=colo.s[5],border=0)

lines(new.dat$Day[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Monoculture")],new.dat$pred[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Monoculture")],lwd=2,lty=2,col=colo.l[2])
points(new.dat$Day[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Monoculture")],log(Eco2.data$Bact_per_ul[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Monoculture")]),lwd=2,col=alpha(colo.l[2],0.4),pch=17)
lines(new.dat$Day[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Resource.alone")],new.dat$pred[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Resource.alone")],lwd=2,lty=2,col=colo.l[1])
points(new.dat$Day[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Resource.alone")],log(Eco2.data$Bact_per_ul[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Resource.alone")]),lwd=2,col=alpha(colo.l[1],0.4),pch=17)
lines(new.dat$Day[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Competition")],new.dat$pred[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Competition")],lwd=2,lt=2,col=colo.l[3])
points(new.dat$Day[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Competition")],log(Eco2.data$Bact_per_ul[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Competition")]),lwd=2,col=alpha(colo.l[3],0.4),pch=17)
lines(new.dat$Day[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Facilitation")],new.dat$pred[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Facilitation")],lwd=2,lty=2,col=colo.l[4])
points(new.dat$Day[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Facilitation")],log(Eco2.data$Bact_per_ul[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Facilitation")]),lwd=2,col=alpha(colo.l[4],0.4),pch=17)
lines(new.dat$Day[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Predation")],new.dat$pred[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Predation")],lwd=2,lty=2,col=colo.l[5])
points(new.dat$Day[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Predation")],log(Eco2.data$Bact_per_ul[which(new.dat$Composition == "Tet" & new.dat$Interaction=="Predation")]),lwd=2,col=alpha(colo.l[5],0.4),pch=17)


#..Statistical models for Tetrahymena density in Ecosystem 2 - (AIC-based simplification procedure - see Methods)
#Linear mixed model 

mod2.3 = with(Eco2.data.tet,lme(log(Tet) ~ Interaction*Day,random = ~ Day|Replicate,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim")))
summary(mod2.3)$tTable
summary(mod2.3)$AIC

mod2.3.1 = with(Eco2.data.tet,update(mod2.3, ~. - Interaction:Day))
summary(mod2.3.1)$AIC
anova(mod2.3,mod2.3.1)


mod2.3.2 = with(Eco2.data.tet,update(mod2.3.1, ~. - Day))
summary(mod2.3.2)$AIC
anova(mod2.3.1,mod2.3.2)

#Rerun final model with REML
mod2.4 = with(Eco2.data.tet,lme(log(Tet) ~ Interaction+Day,random = ~ Day|Replicate,method="REML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim")))
summary(mod2.4)$tTable
coefplot2(mod2.4,vertical=T)
dev.off()
intervals(mod2.4)

#Extract and plot model prediction + CI (http://stackoverflow.com/questions/14358811/extract-prediction-band-from-lme-fit)
X2 = rep(seq(3,14,len=4),times=1,each=40)
new.dat2 = data.frame(Day=X2,Interaction=Eco2.data.tet$Interaction)
new.dat2$pred = predict(mod2.4,newdata=new.dat2,level=0)
Designmat2 = model.matrix(eval(eval(mod2.4$call$fixed)[-2]),new.dat2[-ncol(new.dat2)])
predvar2 = diag(Designmat2 %*% mod2.4$varFix %*% t(Designmat2))
new.dat2$SE = sqrt(predvar2)

#Figure 2 c)

colo.s = c("#99999930", "#E69F0030", "#56B4E930", "#009E7330","#CC79A730")
colo.l = c("#999999", "#E69F00", "#56B4E9", "#009E73","#CC79A7")
#Figure frame
with(Eco2.data.tet,plot(Day,log(Tet),type="n",xaxt="n",yaxt="n"))
axis(side = 2, at =c(-5.6,-3,-0.72,0.096),labels=c(round(exp(-5.6),2),round(exp(-3),2),round(exp(-0.72),2),round(exp(0.096),2))) #Use with dilution 0.33
#axis(side = 2, at =c(0,0.5,1,2),labels=c(round(exp(0),0),round(exp(0.5),0),round(exp(1),0),round(exp(2),0)))
axis(side = 1, at =c(0,3,7,10,14))

#Draw Shadings first and then lines (to avoid having shadings hiding lines)
x11 = new.dat2$Day[which(new.dat2$Interaction=="Monoculture")]
y1.11 = new.dat2$pred[which(new.dat2$Interaction=="Monoculture")]+new.dat2$SE[which(new.dat2$Interaction=="Monoculture")]
y2.11 = new.dat2$pred[which(new.dat2$Interaction=="Monoculture")]-new.dat2$SE[which(new.dat2$Interaction=="Monoculture")]
polygon(c(x11,rev(x11)),c(y2.11,rev(y1.11)),col=colo.s[2],border=0)

x12 = new.dat2$Day[which(new.dat2$Interaction=="Resource.alone")]
y1.12 = new.dat2$pred[which(new.dat2$Interaction=="Resource.alone")]+new.dat2$SE[which(new.dat2$Interaction=="Resource.alone")]
y2.12 = new.dat2$pred[which(new.dat2$Interaction=="Resource.alone")]-new.dat2$SE[which(new.dat2$Interaction=="Resource.alone")]
polygon(c(x12,rev(x12)),c(y2.12,rev(y1.12)),col=colo.s[1],border=0)

x13 = new.dat2$Day[which(new.dat2$Interaction=="Competition")]
y1.13 = new.dat2$pred[which(new.dat2$Interaction=="Competition")]+new.dat2$SE[which(new.dat2$Interaction=="Competition")]
y2.13 = new.dat2$pred[which(new.dat2$Interaction=="Competition")]-new.dat2$SE[which(new.dat2$Interaction=="Competition")]
polygon(c(x13,rev(x13)),c(y2.13,rev(y1.13)),col=colo.s[3],border=0)

x14 = new.dat2$Day[which(new.dat2$Interaction=="Facilitation")]
y1.14 = new.dat2$pred[which(new.dat2$Interaction=="Facilitation")]+new.dat2$SE[which(new.dat2$Interaction=="Facilitation")]
y2.14 = new.dat2$pred[which(new.dat2$Interaction=="Facilitation")]-new.dat2$SE[which(new.dat2$Interaction=="Facilitation")]
polygon(c(x14,rev(x14)),c(y2.14,rev(y1.14)),col=colo.s[4],border=0)

x15 = new.dat2$Day[which(new.dat2$Interaction=="Predation")]
y1.15 = new.dat2$pred[which(new.dat2$Interaction=="Predation")]+new.dat2$SE[which(new.dat2$Interaction=="Predation")]
y2.15 = new.dat2$pred[which(new.dat2$Interaction=="Predation")]-new.dat2$SE[which(new.dat2$Interaction=="Predation")]
polygon(c(x15,rev(x15)),c(y2.15,rev(y1.15)),col=colo.s[5],border=0)

lines(new.dat2$Day[which(new.dat2$Interaction=="Monoculture")],new.dat2$pred[which(new.dat2$Interaction=="Monoculture")],lwd=2,lty=2,col=colo.l[2])
points(Eco2.data.tet$Day[which(new.dat2$Interaction=="Monoculture")],log(Eco2.data.tet$Tet[which(new.dat2$Interaction=="Monoculture")]),lwd=2,col=alpha(colo.l[2],0.4),pch=17)
lines(new.dat2$Day[which(new.dat2$Interaction=="Resource.alone")],new.dat2$pred[which(new.dat2$Interaction=="Resource.alone")],lwd=2,lty=2,col=colo.l[1])
points(Eco2.data.tet$Day[which(new.dat2$Interaction=="Resource.alone")],log(Eco2.data.tet$Tet[which(new.dat2$Interaction=="Resource.alone")]),lwd=2,col=alpha(colo.l[1],0.4),pch=17)
lines(new.dat2$Day[which(new.dat2$Interaction=="Competition")],new.dat2$pred[which(new.dat2$Interaction=="Competition")],lwd=2,lty=2,col=colo.l[3])
points(Eco2.data.tet$Day[which(new.dat2$Interaction=="Competition")],log(Eco2.data.tet$Tet[which(new.dat2$Interaction=="Competition")]),lwd=2,col=alpha(colo.l[3],0.4),pch=17)
lines(new.dat2$Day[which(new.dat2$Interaction=="Facilitation")],new.dat2$pred[which(new.dat2$Interaction=="Facilitation")],lwd=2,lty=2,col=colo.l[4])
points(Eco2.data.tet$Day[which(new.dat2$Interaction=="Facilitation")],log(Eco2.data.tet$Tet[which(new.dat2$Interaction=="Facilitation")]),lwd=2,col=alpha(colo.l[4],0.4),pch=17)
lines(new.dat2$Day[which(new.dat2$Interaction=="Predation")],new.dat2$pred[which(new.dat2$Interaction=="Predation")],lwd=2,lty=2,col=colo.l[5])
points(Eco2.data.tet$Day[which(new.dat2$Interaction=="Predation")],log(Eco2.data.tet$Tet[which(new.dat2$Interaction=="Predation")]),lwd=2,col=alpha(colo.l[5],0.4),pch=17)



#########################################################################################################
#...Statistical Analysis for Ecosystem 1 

#..Identify probability distribution for Bacteria density in Eco 1
summary(Eco1.data$Bact_per_ul)
with(Eco1.data,plot(density(Bact_per_ul)))
with(Eco1.data,plot(density(log(Bact_per_ul))))
with(Eco1.data,qqp(Bact_per_ul,"norm"))
with(Eco1.data,qqp(Bact_per_ul,"lnorm"))
poisson <- with(Eco1.data,fitdistr(Bact_per_ul, "Poisson"))
with(Eco1.data,qqp(Bact_per_ul, "pois",poisson$estimate))
gamma = with(Eco1.data,fitdistr(Bact_per_ul, "gamma")) #only for value larger than zero
with(Eco1.data,qqp(Bact_per_ul, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]]))


#Linear mixed model - (AIC-based simplification procedure - see Methods)
mod3.3 = with(Eco1.data,lme(Bact_per_ul ~ Interaction*Day,random = ~ Day|Replicate,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim")))
summary(mod3.3)$tTable
summary(mod3.3)$AIC

mod3.3.1 = with(Eco1.data,update(mod3.3, ~. - Interaction:Day))
summary(mod3.3.1)$tTable
summary(mod3.3.1)$AIC
anova(mod3.3,mod3.3.1)

#run final model with REML
mod3.4 = with(Eco1.data,lme(Bact_per_ul ~ Interaction*Day,random = ~ Day|Replicate,method="REML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim")))
summary(mod3.4)$tTable
summary(mod3.4)$AIC

#Extract and plot model prediction + CI (http://stackoverflow.com/questions/14358811/extract-prediction-band-from-lme-fit)
X3 = rep(seq(3,14,len=4),times=1,each=40)
new.dat3 = data.frame(Day=X3,Interaction=Eco1.data$Interaction)
new.dat3$pred = predict(mod3.4,newdata=new.dat3,level=0)
Designmat3 = model.matrix(eval(eval(mod3.4$call$fixed)[-2]),new.dat3[-ncol(new.dat3)])
predvar3 = diag(Designmat3 %*% mod3.4$varFix %*% t(Designmat3))
new.dat3$SE = sqrt(predvar3)


#Figure 2 a)

colo.s = c("#99999930", "#E69F0030", "#56B4E930", "#009E7330","#CC79A730")
colo.l = c("#999999", "#E69F00", "#56B4E9", "#009E73","#CC79A7")
#Figure frame
with(Eco1.data,plot(Day,Bact_per_ul,type="n",xaxt="n",yaxt="n",log="y"))
axis(side = 2, at =c(1400,7000,40000,120000,230000)) # to use with Dil 0.33
#axis(side = 2, at =c(7000,40000,120000,200000,360000))
axis(side = 1, at =c(0,3,7,10,14))#   ,labels=ylabels[[i]])

#Draw Shadings first and then lines (to avoid having shadings hiding lines)
x12 = new.dat3$Day[which(new.dat3$Interaction=="Monoculture")]
y1.12 = new.dat3$pred[which(new.dat3$Interaction=="Monoculture")]+new.dat3$SE[which(new.dat3$Interaction=="Monoculture")]
y2.12 = new.dat3$pred[which(new.dat3$Interaction=="Monoculture")]-new.dat3$SE[which(new.dat3$Interaction=="Monoculture")]
polygon(c(x12,rev(x12)),c(y2.12,rev(y1.12)),col=colo.s[2],border=0)

x13 = new.dat3$Day[which(new.dat3$Interaction=="Resource.alone")]
y1.13 = new.dat3$pred[which(new.dat3$Interaction=="Resource.alone")]+new.dat3$SE[which(new.dat3$Interaction=="Resource.alone")]
y2.13 = new.dat3$pred[which(new.dat3$Interaction=="Resource.alone")]-new.dat3$SE[which(new.dat3$Interaction=="Resource.alone")]
polygon(c(x13,rev(x13)),c(y2.13,rev(y1.13)),col=colo.s[1],border=0)

x14 = new.dat3$Day[which(new.dat3$Interaction=="Competition")]
y1.14 = new.dat3$pred[which(new.dat3$Interaction=="Competition")]+new.dat3$SE[which(new.dat3$Interaction=="Competition")]
y2.14 = new.dat3$pred[which(new.dat3$Interaction=="Competition")]-new.dat3$SE[which(new.dat3$Interaction=="Competition")]
polygon(c(x14,rev(x14)),c(y2.14,rev(y1.14)),col=colo.s[3],border=0)

x15 = new.dat3$Day[which(new.dat3$Interaction=="Facilitation")]
y1.15 = new.dat3$pred[which(new.dat3$Interaction=="Facilitation")]+new.dat3$SE[which(new.dat3$Interaction=="Facilitation")]
y2.15 = new.dat3$pred[which(new.dat3$Interaction=="Facilitation")]-new.dat3$SE[which(new.dat3$Interaction=="Facilitation")]
polygon(c(x15,rev(x15)),c(y2.15,rev(y1.15)),col=colo.s[4],border=0)

x16 = new.dat3$Day[which(new.dat3$Interaction=="Predation")]
y1.16 = new.dat3$pred[which(new.dat3$Interaction=="Predation")]+new.dat3$SE[which(new.dat3$Interaction=="Predation")]
y2.16 = new.dat3$pred[which(new.dat3$Interaction=="Predation")]-new.dat3$SE[which(new.dat3$Interaction=="Predation")]
polygon(c(x16,rev(x16)),c(y2.16,rev(y1.16)),col=colo.s[5],border=0)

lines(new.dat3$Day[which(new.dat3$Interaction=="Monoculture")],new.dat3$pred[which(new.dat3$Interaction=="Monoculture")],lwd=2,col=colo.l[2])
points(Eco1.data$Day[which(new.dat3$Interaction=="Monoculture")],Eco1.data$Bact_per_ul[which(new.dat3$Interaction=="Monoculture")],lwd=2,col=alpha(colo.l[2],0.4),pch=16)
lines(new.dat3$Day[which(new.dat3$Interaction=="Resource.alone")],new.dat3$pred[which(new.dat3$Interaction=="Resource.alone")],lwd=2,col=colo.l[1])
points(Eco1.data$Day[which(new.dat3$Interaction=="Resource.alone")],Eco1.data$Bact_per_ul[which(new.dat3$Interaction=="Resource.alone")],lwd=2,col=alpha(colo.l[1],0.4),pch=16)
lines(new.dat3$Day[which(new.dat3$Interaction=="Competition")],new.dat3$pred[which(new.dat3$Interaction=="Competition")],lwd=2,col=colo.l[3])
points(Eco1.data$Day[which(new.dat3$Interaction=="Competition")],Eco1.data$Bact_per_ul[which(new.dat3$Interaction=="Competition")],lwd=2,col=alpha(colo.l[3],0.4),pch=16)
lines(new.dat3$Day[which(new.dat3$Interaction=="Facilitation")],new.dat3$pred[which(new.dat3$Interaction=="Facilitation")],lwd=2,col=colo.l[4])
points(Eco1.data$Day[which(new.dat3$Interaction=="Facilitation")],Eco1.data$Bact_per_ul[which(new.dat3$Interaction=="Facilitation")],lwd=2,col=alpha(colo.l[4],0.4),pch=16)
lines(new.dat3$Day[which(new.dat3$Interaction=="Predation")],new.dat3$pred[which(new.dat3$Interaction=="Predation")],lwd=2,col=colo.l[5])
points(Eco1.data$Day[which(new.dat3$Interaction=="Predation")],Eco1.data$Bact_per_ul[which(new.dat3$Interaction=="Predation")],lwd=2,col=alpha(colo.l[5],0.4),pch=16)




#########################################################################################################
#...Exploratory figures

#Figure 3

#a) Downstream bacteria density as a function of upstream bacteria density

pdf(paste0(graphpath,"Figure3_ver3_(none-diluted).pdf"),width=8,height=8)

#determine min and max for y axes (change as a function of whether you are ploting diluted or none-diluted systems)
min = 0
max = 0
if(Eco1.data$Dilution[1]==0) {
  min=20000
  max= 500000} else{min=10000
  max=200000}

plot(Eco1.data$Bact_per_ul[which(Eco1.data$Structure=="Consumer")],Eco2.data$Bact_per_ul[which(Eco2.data$Structure=="Consumer")],ylim=c(min,max),pch=16,xlab="Upstream bacteria density (ind/uL)", ylab="Downstream bacteria density (ind/uL)",xaxt="n",yaxt="n",type="n")
if(Eco1.data$Dilution[1]==0) {
  axis(side = 2, at =c(20000,100000,200000,300000,400000,500000),labels=c("20000","100000","200000","300000","400000","500000")) 
  axis(side = 1, at =c(0,50000,100000,150000,200000,250000,300000,350000))} else{
    axis(side = 2, at =c(10000,50000,100000,150000,200000),labels=c("10000","50000","100000","150000","200000")) 
    axis(side = 1, at =c(0,50000,100000,150000,200000))
  }
points(Eco1.data$Bact_per_ul[which(Eco1.data$Structure=="Resource")],Eco2.data$Bact_per_ul[which(Eco2.data$Structure=="Resource")],pch=16,col=alpha("blue",0.8))
abline(lm(Eco2.data$Bact_per_ul[which(Eco2.data$Structure=="Resource")]~Eco1.data$Bact_per_ul[which(Eco1.data$Structure=="Resource")]),col="blue",lwd=2)
points(Eco1.data$Bact_per_ul[which(Eco1.data$Structure=="Consumer")],Eco2.data$Bact_per_ul[which(Eco2.data$Structure=="Consumer")],pch=16,col=alpha("red",0.8))
abline(lm(Eco2.data$Bact_per_ul[which(Eco2.data$Structure=="Consumer")]~Eco1.data$Bact_per_ul[which(Eco1.data$Structure=="Consumer")]),col="red",lwd=2)

#a) Downstream Tetrahymena density as a function of upstream bacteria density

#determine min and max for y axes (change as a function of whether you are ploting diluted or none-diluted systems)
min = 0
max = 0
if(Eco1.data$Dilution[1]==0) {
  min=0
  max= 6} else{min=0
  max=1.1}

plot(Eco1.data$Bact_per_ul[which(Eco1.data$Structure=="Consumer")],Eco2.data$Tet[which(Eco2.data$Structure=="Consumer")],pch=16,ylim=c(min,max),xlab="Upstream bacteria density (ind/uL)", ylab="Downstream Tetrahymena density (ind/uL)",col=alpha("red",0.8))
abline(lm(Eco2.data$Tet[which(Eco2.data$Structure=="Consumer")]~Eco1.data$Bact_per_ul[which(Eco1.data$Structure=="Consumer")]),col="red",lwd=2)


dev.off()

#.........Ecosystem 1 protist dynamics

#Resource alone
lineplot.CI(Day,Bact_per_ul,col="lightgreen",xlab="",ylab="Ind./ul",data= Eco1.data[Eco1.data$Dilution==0 & Eco1.data$Interaction=="Resource.alone",])
#Monoculture
lineplot.CI(Day,Col,col="lightblue",xlab="",ylab="",data= Eco1.data[Eco1.data$Dilution==0 & Eco1.data$Interaction=="Monoculture",])
#Competition
lineplot.CI(Day,Col,col="lightblue",xlab="",ylab="",ylim=c(0,0.52),data= Eco1.data[Eco1.data$Dilution==0 & Eco1.data$Interaction=="Competition",])
par(new=T)
lineplot.CI(Day,Pau,col="pink",xlab="",ylab="",ylim=c(0,0.52),data= Eco1.data[Eco1.data$Dilution==0 & Eco1.data$Interaction=="Competition",])
#Facilitation
lineplot.CI(Day,Col,col="lightblue",xlab="",ylab="",ylim=c(0,1.7),data= Eco1.data[Eco1.data$Dilution==0 & Eco1.data$Interaction=="Facilitation",])
par(new=T)
lineplot.CI(Day,Eug,col="darkgreen",xlab="",ylim=c(0,1.7),ylab="",data= Eco1.data[Eco1.data$Dilution==0 & Eco1.data$Interaction=="Facilitation",])
#Predation
lineplot.CI(Day,Col,col="lightblue",xlab="",ylim=c(0,3),ylab="",data= Eco1.data[Eco1.data$Dilution==0 & Eco1.data$Interaction=="Predation",])
par(new=T)
lineplot.CI(Day,Daphnia,col="red",xlab="",ylim=c(0,3),ylab="",data= Eco1.data[Eco1.data$Dilution==0 & Eco1.data$Interaction=="Predation",])

#.........Ecosystem 1 Colpidium density for each treatment
colo.l = c("#999999", "#E69F00", "#56B4E9", "#009E73","#CC79A7")
lineplot.CI(Day,Col,col=colo.l[2],xlab="",ylim=c(0,1.7),ylab="Colpidium ind/ul upstream ecosystems",data= Eco1.data[Eco1.data$Interaction=="Monoculture",])
par(new=T)
lineplot.CI(Day,Col,col=colo.l[3],xlab="",ylab="",ylim=c(0,1.7),data= Eco1.data[Eco1.data$Interaction=="Competition",])
par(new=T)
lineplot.CI(Day,Col,col=colo.l[4],xlab="",ylab="",ylim=c(0,1.7),data= Eco1.data[Eco1.data$Interaction=="Facilitation",])
par(new=T)
lineplot.CI(Day,Col,col=colo.l[5],xlab="",ylim=c(0,1.7),ylab="",data= Eco1.data[Eco1.data$Interaction=="Predation",])


#..........Ecosystem 1 - Bacteria FSC.A over time * Community structure


# Along the X-axis is the FSC(Forward SCatter) parameter. 
# This parameter is a measurement of the amount of the laser beam that passes around the cell. 
# This gives a relative size for the cell. 
# If one uses a known control or standard such as beads with a known size, 
# one can determine the relative size of the cells based on the size of the control or standard. 
# Along the Y-axis is the SSC(Side SCatter) parameter. 
# This parameter is a measurement of the amount of the laser beam that bounces off of particulates 
# inside of the cell. Any laser light of the same wavelength of the laser(488nm) 
# that bounces off of things inside the cell is collected and recorded as SSC.
# Using a bead of known size one can determine the size of a population based on the FSC parameters

#Dilution 0
lineplot.CI(Day,Mean.FSC.A,group=Interaction,log="y",xlab="",ylab="Bacteria FSC.A",main="Ecosystem 1 - Bacteria FSC.A - Dilution 0",data=Eco1.data[Eco1.data$Dilution==0 & Eco1.data$Interaction!="Predation" & Eco1.data$Interaction!="Facilitation",])

#Effect of Tetrahymena presence on downstream bacteria coefficient of variation 
#With Tetrahymena 
(sd(Eco2.data$Bact_per_ul[which(Eco2.data$Structure=="Consumer")])/mean(Eco2.data$Bact_per_ul[which(Eco2.data$Structure=="Consumer")]))*100
#Without Tetrahymena 
(sd(Eco2.data$Bact_per_ul[which(Eco2.data$Structure=="Resource")])/mean(Eco2.data$Bact_per_ul[which(Eco2.data$Structure=="Resource")]))*100
