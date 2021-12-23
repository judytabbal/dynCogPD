# This script performs the behavioral analyses of the Simon task data
# A first part consists in overall behavioral analyses (averaged RT and accuracy)
# The second part performs distributional analyses for conditional accuracy functions and delta plots
# This script needs the distrib function to run properly

# duprez.joan@gmail.com


library(ggplot2)
library(xlsx)
library(tidyverse)
library(plotrix)
library(gridExtra)
library(Hmisc)
library(lme4)
library(car)
library(MuMIn)
library(redres)
library(afex)
library(ggsignif)
library(DescTools)


data = read.csv("behavior.csv", sep = "")
colnames(data)=c("num", "DIN", "condition", "RT","accuracy", "subname", "group", "NA")

# Create a n column

subnames = unique(data$subname)
subn = data.frame(
  subname = subnames, 
  ntot = seq(from = 1, to = length(subnames), by = 1))

data <- data %>% right_join(subn, by=c("subname"))

# Remove excluded participants

data = droplevels(data[which(data$subname != '04GC_pilot' & data$subname != '03LM'),])

# Overall data cleaning
ndata = length(data$RT) # used to compute % of removed data
data<-droplevels(data[which(data$RT>200 & data$RT<1500),]) # removes trials with RT<100 & >1500
thresh<-aggregate(data$RT,list(data$ntot),function(x) mean(x)+3*sd(x)) # creates a RT threshold of mean + 3sd
data<-merge(data,thresh,by.x="ntot",by.y="Group.1", all=TRUE) # merges the thresh column with data
data<-droplevels(data[data$RT<data$x,]) # removes trials above threshold
pct_outliers<-(((ndata-length(data$RT))/(ndata)))*100 # Computes the % of removed trials
dataBR=droplevels(data[which(data$accuracy!="0"),])

####                                    ####
####                                    ####
####  Overall RT and accuracy analyses  ####
####                                    ####                              
####                                    ####

# Average RT correct trials only
mean_RTn<-with(dataBR,aggregate(RT,list(dataBR$group,dataBR$ntot,dataBR$condition),mean)) # average RT per participant and condition
mean_RT<-with(dataBR,aggregate(RT,list(dataBR$group,dataBR$condition),mean)) #  average RT per group and condition
mean_RT<-cbind.data.frame(mean_RT,sdRT=with(mean_RTn,aggregate(x,list(Group.1, Group.3),sd))$x)  #  standard deviation
mean_RT<-cbind.data.frame(mean_RT,seRT=with(mean_RTn,aggregate(x,list(Group.1, Group.3),std.error))$x) #standard error (needs plotrix)   
colnames(mean_RTn)<-c("group","n", "condition", "RT")     
colnames(mean_RT)<-c("group","condition","RT","sdRT","seRT")
mean_RT$condition[which(mean_RT$condition=="INCongruent")]<-c("Incongruent") # rename

# Same with accuracy
mean_BRn<-with(data,aggregate(accuracy,list(data$group,data$ntot,data$condition),mean))
mean_BR<-with(data,aggregate(accuracy,list(data$group,data$condition),mean)) 
mean_BR<-cbind.data.frame(mean_BR,sdAcc=with(mean_BRn,aggregate(x,list(Group.1, Group.3),sd))$x)
mean_BR<-cbind.data.frame(mean_BR,seAcc=with(mean_BRn,aggregate(x,list(Group.1, Group.3),std.error))$x) 
colnames(mean_BRn)<-c("group","n", "condition", "Accuracy")     
colnames(mean_BR)<-c("group","condition","Accuracy","sd_Accuracy","se_Accuracy")
mean_BR$condition[which(mean_BR$condition=="1")]<-c("Congruent")
mean_BR$condition[which(mean_BR$condition=="INCongruent")]<-c("Incongruent")



# Plot
g.left<-ggplot(mean_RT, aes(x=condition, y=RT, group=group))+
  geom_line(aes(color=group), size = 0.5) +  
  scale_color_manual(values=c("#1f77b4", "#ff7f0e"))+
  geom_errorbar(aes(ymin=RT-seRT, ymax=RT+seRT), colour="black", width=0.1)+
  geom_point(aes(shape=group, colour=group), size=4)+  
  ggtitle("A")+
  theme(legend.justification=c(1,0),legend.position=c(0.25,0.70),
        plot.title = element_text(face = "bold"),
        axis.title.y=element_text(vjust=1.5),
        axis.title.x = element_blank(),
        text=element_text(colour="black", size=15),
        panel.background = element_rect(fill=NA, colour='black', size=1.2))+
  scale_y_continuous(name="RT (ms)")


g.right<-ggplot(mean_BR, aes(x=condition, y=Accuracy, group=group))+
  geom_line(aes(color=group), size = 0.5) +  
  scale_color_manual(values=c("#1f77b4", "#ff7f0e"))+
  geom_errorbar(aes(ymin=Accuracy-se_Accuracy, ymax=Accuracy+se_Accuracy), colour="black", width=0.1)+
  geom_point(aes(shape=group, colour=group), size=4)+   
  ggtitle("B")+
  theme(legend.justification=c(1,0),legend.position="none",
        plot.title = element_text(face = "bold"),
        axis.title.y=element_text(vjust=1.5),
        axis.title.x = element_blank(),
        text=element_text(colour="black", size=15),
        panel.background = element_rect(fill=NA, colour='black', size=1.2))+
  scale_y_continuous(name="Accuracy")

print(grid.arrange(g.left,g.right, heights=c(1/2,1/2)))


## Stats congruence and group effect

# RT
dataBR$condition = as.factor(dataBR$condition)
dataBR$group = as.factor(dataBR$group)

RTmod = lmer(log(RT)~condition*group + (condition|subname), data=dataBR)
Anova(RTmod)
r.squaredGLMM(RTmod)
launch_redres(RTmod)

# Accuracy
data$condition = as.factor(data$condition)
data$group = as.factor(data$group)
Accmod = glmer(accuracy~condition*group + (condition|subname), family=binomial,data=data, glmerControl(optimizer="bobyqa"))
Anova(Accmod)
r.squaredGLMM(Accmod)


####                                                  ####
####                                                  ####
####  Conditional accuracy functions and delta plots  ####
####                                                  ####                              
####                                                  #### make sure the distrib function is in the path



## CAF

nquantiles=7

## loop through subjects
subs = unique(data$ntot)

r<-1
matC<-matrix(NA,nquantiles*length(subs),6) ## Congruent
matNC<-matrix(NA,nquantiles*length(subs),6) ## Inongruent


for(i in 1:length(subs)){
  nC<-distrib(data, subs[i],"Congruent", nquantiles)
  matC[r:(r+nquantiles-1),]<-cbind(round(nC,3), rep(unique(data[which(data$ntot == subs[i]),]$ntot), nquantiles),
                                   rep(unique(data[which(data$ntot == subs[i]),]$ntot), nquantiles), 
                                   rep(unique(data[which(data$ntot == subs[i]),]$subname), nquantiles),
                                   rep(unique(data[which(data$ntot == subs[i]),]$group), nquantiles))
  
  nNC<-distrib(data, subs[i],"INCongruent", nquantiles)
  matNC[r:(r+nquantiles-1),]<-cbind(round(nNC,3), rep(unique(data[which(data$ntot == subs[i]),]$ntot), nquantiles),
                                    rep(unique(data[which(data$ntot == subs[i]),]$ntot), nquantiles), 
                                    rep(unique(data[which(data$ntot == subs[i]),]$subname), nquantiles),
                                    rep(unique(data[which(data$ntot == subs[i]),]$group), nquantiles))
  r<-r+nquantiles
}

C<-data.frame(cbind(n=sort(rep(c(1:length(subs)),nquantiles)), 
                    condition= "Congruent", quantile=1:nquantiles), matC)
NC<-data.frame(cbind(n=sort(rep(c(1:length(subs)),nquantiles)), 
                     condition= "Incongruent", quantile=1:nquantiles), matNC)

TotalBRTR<-rbind(C,NC) # final table with all trials

names(TotalBRTR)=c("ntot","condition","quantile","RT","accuracy", "n", "ntot2", "subname", "group")
TotalBRTR$RT=as.numeric(TotalBRTR$RT)
TotalBRTR$accuracy=as.numeric(TotalBRTR$accuracy)

## Average per bin

mTotalBRTR<-aggregate(list(TotalBRTR$RT,TotalBRTR$accuracy),
                      list(TotalBRTR$quantile,TotalBRTR$condition, TotalBRTR$group),
                      function(x) round(mean(x),3)) # seuil RT pour chaque n
sdTotalBRTR<-aggregate(list(TotalBRTR$RT,TotalBRTR$accuracy),
                       list(TotalBRTR$quantile,TotalBRTR$condition, TotalBRTR$group),
                       function(x) round(sd(x),3)) # seuil RT pour chaque n
seTotalBRTR<-aggregate(list(TotalBRTR$RT,TotalBRTR$accuracy),
                       list(TotalBRTR$quantile,TotalBRTR$condition, TotalBRTR$group),
                       function(x) round(sd(x)/sqrt(max(data$ntot)),3)) # seuil RT pour chaque n
names(sdTotalBRTR)=c("quantile","condition","group","RT", "accuracy")
names(seTotalBRTR)=c("quantile","condition","group","RT", "accuracy")
mTotalBRTR<-cbind(mTotalBRTR,sdTotalBRTR$RT,sdTotalBRTR$accuracy,seTotalBRTR$RT,seTotalBRTR$accuracy) ## Tableau final moyenne par condition, gain et quantile
names(mTotalBRTR)=c("quantile","condition","group","RT","accuracy","sdRT","sdBR","seRT","seBR")

## Plot CAF

ggplot(mTotalBRTR[which(mTotalBRTR$condition=="Incongruent"),], aes(x=RT, y=accuracy, group=group))+
  coord_cartesian(ylim = c(0.5, 1), xlim=c(400,850))+
  scale_color_manual(values=c("#1f77b4", "#ff7f0e"))+
  geom_errorbar(aes(ymin=accuracy-seBR, ymax=accuracy+seBR),  width=6)+
  geom_line(aes(color=group), size = 1) +
  geom_errorbarh(aes(xmin=RT-seRT, xmax=RT+seRT), height=0.005)+  
  geom_point(aes(shape=group, color=group), size=4,  fill="white")  +    
  scale_x_continuous(name="RT (ms)",breaks = seq(400, 850, 100))+
  theme(axis.title.y=element_text(vjust=1.5),axis.title.x=element_text(vjust=0.1),
        legend.justification=c(1,0),legend.position=c(0.99,0.01), 
        text=element_text(colour="black", size=20),
        panel.background = element_rect(fill=NA, colour='black', size=1.2))+
  scale_y_continuous(name="Incongruent accuracy")

ggplot(mTotalBRTR[which(mTotalBRTR$condition=="Congruent"),], aes(x=RT, y=accuracy, group=group))+
  coord_cartesian(ylim = c(0.5, 1), xlim=c(400,850))+
  geom_line(aes(color=group), size = 1) +
  scale_color_manual(values=c("#1f77b4", "#ff7f0e"))+
  geom_errorbar(aes(ymin=accuracy-seBR, ymax=accuracy+seBR),  width=6)+
  geom_errorbarh(aes(xmin=RT-seRT, xmax=RT+seRT), height=0.005)+  
  geom_point(aes(shape=group, color=group), size=4,  fill="white")  +    
  scale_x_continuous(name="RT (ms)",breaks = seq(400, 850, 100))+
  theme(axis.title.y=element_text(vjust=1.5),axis.title.x=element_text(vjust=0.1),
        legend.justification=c(1,0),legend.position=c(0.99,0.01), 
        text=element_text(colour="black", size=20),
        panel.background = element_rect(fill=NA, colour='black', size=1.2))+
  scale_y_continuous(name="Congruent accuracy")

# Test violinplots
dodge <- position_dodge(width = 0.5)
caf=ggplot(TotalBRTR[which(TotalBRTR$quantile == 1 & TotalBRTR$condition == "Incongruent"),], 
           aes(x=group, y=accuracy, fill=group)) + 
  geom_violin(width=0.75, position = dodge, alpha = .5)+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  scale_fill_manual(values=c("#1f77b4", "#ff7f0e"))+
  scale_y_continuous(name="First bin incongruent accuracy")+
  ggtitle("C")+
  theme(axis.title.y=element_text(vjust=1.5),axis.title.x=element_text(vjust=0.1),
        legend.justification=c(1,0),legend.position="none", 
        plot.title = element_text(face = "bold"),
        text=element_text(colour="black", size=15), 
        panel.background = element_rect(fill=NA, colour='black', size=1.2))

# Stat group effect

t.test(TotalBRTR[which(TotalBRTR$quantile == 1 & TotalBRTR$condition == "Incongruent" & TotalBRTR$group == "HC"),]$accuracy,
       TotalBRTR[which(TotalBRTR$quantile == 1 & TotalBRTR$condition == "Incongruent" & TotalBRTR$group == "PD"),]$accuracy)

CohenD(TotalBRTR[which(TotalBRTR$quantile == 1 & TotalBRTR$condition == "Incongruent" & TotalBRTR$group == "HC"),]$accuracy,
       TotalBRTR[which(TotalBRTR$quantile == 1 & TotalBRTR$condition == "Incongruent" & TotalBRTR$group == "PD"),]$accuracy, pooled=TRUE)


## Delta plots

## Boucle appliquant la fonction ? tous les participants
r<-1
matC<-matrix(NA,nquantiles*n,4) ## Congruent
matNC<-matrix(NA,nquantiles*n,4) ## Non Congruent

subs = unique(data$ntot)

for(i in 1:n){
  nC<-distrib_delta(dataBR, subs[i],"Congruent", nquantiles)
  matC[r:(r+nquantiles-1),]<-cbind(round(nC,3),
                                   rep(unique(dataBR[which(dataBR$ntot == subs[i]),]$ntot), nquantiles), 
                                   rep(unique(dataBR[which(dataBR$ntot == subs[i]),]$subname), nquantiles),
                                   rep(unique(dataBR[which(dataBR$ntot == subs[i]),]$group), nquantiles))
  
  nNC<-distrib_delta(dataBR, subs[i],"INCongruent", nquantiles)
  matNC[r:(r+nquantiles-1),]<-cbind(round(nNC,3),
                                    rep(unique(dataBR[which(dataBR$ntot == subs[i]),]$ntot), nquantiles), 
                                    rep(unique(dataBR[which(dataBR$ntot == subs[i]),]$subname), nquantiles),
                                    rep(unique(dataBR[which(dataBR$ntot == subs[i]),]$group), nquantiles))
  r<-r+nquantiles
}


## Regroupement des données par congruence 
C<-data.frame(cbind(congruence=1, quantile=1:nquantiles), matC)
names(C)=c("congruence","quantile","RT", "n", "subname","group")
C$RT=as.numeric(C$RT)
# C = droplevels(C[which(C$subname != '03CC_pilot'),])
# C = droplevels(C[which(C$subname != '03LM'),])


NC<-data.frame(cbind(congruence=0, quantile=1:nquantiles), matNC)
names(NC)=c("congruence","quantile","RT", "n", "subname","group")
NC$RT=as.numeric(NC$RT)
# NC = droplevels(NC[which(NC$subname != '03CC_pilot'),])
# NC = droplevels(NC[which(NC$subname != '03LM'),])

Total<-rbind(C,NC) # Tableau final tous essais
names(Total)=c("congruence","quantile","RT", "n", "subname","group")

## Calcul des TR moyens par points
nTotal<-aggregate(Total$RT,list(n=Total$n,group=Total$group,quantile=Total$quantile),function(x) round(mean(x),3)) # TR du quantile par participant
mTotal<-aggregate(Total$RT,list(Total$quantile, Total$group),function(x) round(mean(x),3)) # TR moyen par quantile
sdTotal<-aggregate(Total$RT,list(Total$quantile, Total$group),function(x) round(sd(x),3)) # seuil TR pour chaque n
seTotal<-aggregate(Total$RT,list(Total$quantile,Total$group),function(x) round(sd(x)/sqrt(max(data$ntot)),3)) # seuil TR pour chaque n
names(nTotal)<-c("n","group","quantile","RT")
names(mTotal)=c("quantile", "group","RT")
names(sdTotal)=c("points", "group","RT")
names(seTotal)=c("points", "group","RT")

## Calcul des Simon moyens par points

datanSimon<-data.frame(n=C$n,quantile=C$quantile,Simon=NC$RT-C$RT, group=C$group)# simon moyen par point et par participant, colonne 9 est le groupe dans les tableaux C et NC


mSimon<-aggregate(datanSimon$Simon,list(datanSimon$quantile, datanSimon$group),function(x) round(mean(x),3)) # Simon moyen par point
sdSimon<-aggregate(datanSimon$Simon,list(datanSimon$quantile, datanSimon$group),function(x) round(sd(x),3)) # Ecart type
seSimon<-aggregate(datanSimon$Simon,list(datanSimon$quantile,datanSimon$group),function(x)round(sd(x)/sqrt(max(data$ntot)),3)) # Erreur standard
names(mSimon)=c("quantile", "group", "delta_RT")
names(sdSimon)=c("quantile","group", "sd_RT")
names(seSimon)=c("quantile","group","se_RT")

mSimon2<-cbind(mSimon,seSimon$se_RT,sdSimon$sd_RT,mTotal$RT,seTotal$RT,sdTotal$RT)
colnames(mSimon2)<-c("quantile", "group", "Delta_RT","se_Delta","sd_Delta","RT","se_RT","sd_RT")


## Plot delta plots

ggplot(mSimon2, aes(x=RT, y=Delta_RT, group=group))+coord_cartesian(xlim=c(400,850))+
  geom_line(aes(color=group), size = 1) +
  scale_color_manual(values=c("#1f77b4", "#ff7f0e"))+
  geom_errorbar(aes(ymin=Delta_RT-se_Delta, ymax=Delta_RT+se_Delta),  width=6)+
  geom_errorbarh(aes(xmin=RT-se_RT, xmax=RT+se_RT), height=0.005)+  
  geom_point(aes(shape=group, color=group), size=4,  fill="white")  +    
  scale_x_continuous(name="RT (ms)",breaks = seq(400, 850, 100))+
  theme(axis.title.y=element_text(vjust=1.5),axis.title.x=element_text(vjust=0.1),legend.justification=c(1,0),legend.position=c(0.99,0.01), text=element_text(colour="black", size=20),panel.background = element_rect(fill=NA, colour='black', size=1.2))+
  scale_y_continuous(name="Congruence effect (ms)")


## Slope calculation
#gp1
# Y
y7pts<-droplevels(datanSimon[datanSimon$quantile==(nquantiles),]) ## Donnees 7eme point 
y6pts<-droplevels(datanSimon[datanSimon$quantile==(nquantiles-1),]) ## Donnees 6eme point
y5pts<-droplevels(datanSimon[datanSimon$quantile==(nquantiles-2),]) ## Donnees 5eme point
y4pts<-droplevels(datanSimon[datanSimon$quantile==(nquantiles-3),]) ## Donnees 4eme point
y3pts<-droplevels(datanSimon[datanSimon$quantile==(nquantiles-4),]) ## Donnees 3eme point
y2pts<-droplevels(datanSimon[datanSimon$quantile==(nquantiles-5),]) ## Donnees 2nd point
y1pts<-droplevels(datanSimon[datanSimon$quantile==(nquantiles-6),]) ## Donnees 1er point

# slopes
y67<-data.frame(n=y7pts$n,npente=6,ypente=y6pts$Simon-y7pts$Simon, group=y7pts$group)
y56<-data.frame(n=y6pts$n,npente=5,ypente=y5pts$Simon-y6pts$Simon, group=y7pts$group)
y45<-data.frame(n=y5pts$n,npente=4,ypente=y4pts$Simon-y5pts$Simon, group=y7pts$group)
y34<-data.frame(n=y4pts$n,npente=3,ypente=y3pts$Simon-y4pts$Simon, group=y7pts$group)
y23<-data.frame(n=y3pts$n,npente=2,ypente=y2pts$Simon-y3pts$Simon, group=y7pts$group)
y12<-data.frame(n=y2pts$n,npente=1,ypente=y1pts$Simon-y2pts$Simon, group=y7pts$group)


y_gp1<-rbind(y12,y23,y34,y45,y56,y67)

# X
x7pts<-droplevels(nTotal[nTotal$quantile==(nquantiles),]) ## Données 7eme point 
x6pts<-droplevels(nTotal[nTotal$quantile==(nquantiles-1),]) ## Données 6eme point
x5pts<-droplevels(nTotal[nTotal$quantile==(nquantiles-2),]) ## Données 5eme point
x4pts<-droplevels(nTotal[nTotal$quantile==(nquantiles-3),]) ## Données 4eme point
x3pts<-droplevels(nTotal[nTotal$quantile==(nquantiles-4),]) ## Données 3eme point
x2pts<-droplevels(nTotal[nTotal$quantile==(nquantiles-5),]) ## Données 2nd point
x1pts<-droplevels(nTotal[nTotal$quantile==(nquantiles-6),]) ## Données 1er point

# slopes

x67<-data.frame(n=x5pts$n,npente=6,xpente=x6pts$RT-x7pts$RT, group=x7pts$group)
x56<-data.frame(n=x5pts$n,npente=5,xpente=x5pts$RT-x6pts$RT, group=x7pts$group)
x45<-data.frame(n=x5pts$n,npente=4,xpente=x4pts$RT-x5pts$RT, group=x7pts$group)
x34<-data.frame(n=x5pts$n,npente=3,xpente=x3pts$RT-x4pts$RT, group=x7pts$group)
x23<-data.frame(n=x5pts$n,npente=2,xpente=x2pts$RT-x3pts$RT, group=x7pts$group)
x12<-data.frame(n=x5pts$n,npente=1,xpente=x1pts$RT-x2pts$RT, group=x7pts$group)

x_gp1<-rbind(x12,x23,x34,x45,x56,x67)


pentes<-data.frame(n=x_gp1$n,npente=x_gp1$npente,pente=round((y_gp1$ypente/x_gp1$xpente),3),group=x_gp1$group)

# Plot last slope

dodge <- position_dodge(width = 0.5)
p= ggplot(pentes[which(pentes$npente == 6),], aes(x=group, y=pente, fill=group)) + 
  geom_violin(width=0.75, position = dodge, alpha = .5)+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  scale_fill_manual(values=c("#1f77b4", "#ff7f0e"))+
  scale_y_continuous(name="Last slope value", limits = c(-1, 0.125))+
  ggtitle("D")+
  geom_signif(comparisons = list(c("HC", "PD")), 
              annotation = c("* (p = 0.023)"), map_signif_level=TRUE)+
  theme(axis.title.y=element_text(vjust=1.5),axis.title.x=element_text(vjust=0.1),
        plot.title = element_text(face = "bold"),
        legend.justification=c(1,0),legend.position="none", text=element_text(colour="black", size=15),
        panel.background = element_rect(fill=NA, colour='black', size=1.2))

# Stat group effect
t.test(list(pentes[which(pentes$group=='HC' & pentes$npente==6),]$pente, 
                 pentes[which(pentes$group=='PD' & pentes$npente==6),]$pente))

CohenD(pentes[which(pentes$group=='HC' & pentes$npente==6),]$pente, pentes[which(pentes$group=='PD' & pentes$npente==6),]$pente, pooled=TRUE)
