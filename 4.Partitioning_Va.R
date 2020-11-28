library(cowplot)
library(ggplot2)
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat) 
}

dir <- "~/Dropbox/dropbox_work/data/GWAS_m_f_fitness/"

##############################
## Prep autosome/X sites for REML
##############################

## Make list1 and list2 files; these are simply a column of predictors situated on autosomes vs X chromosome
lhm.bim <- read.table(paste0(dir,"f3c.lhm.snp.bim"))
list1 <- subset(lhm.bim,V1!=5)$V2  
list2 <- subset(lhm.bim,V1==5)$V2
write.table(list1,paste0(dir,"autosome_X_partition/list1"),quote=F,row.names=F,col.names=F)
write.table(list2,paste0(dir,"autosome_X_partition/list2"),quote=F,row.names=F,col.names=F)

##############################
## Circular genotypic permutations
##############################

## 1,000 circular permutations of list1 and list2 files
set.seed(123)
start.point <- vector("numeric")
for (i in 1:1000){
  chrom.circular <- rep(lhm.bim$V1,2)
  start.point[i] <- sample(1:765764,1)
  lhm.bim$V1_perm <- chrom.circular[start.point[i]:(start.point[i]+765763)]
  list1.perm <- subset(lhm.bim,V1_perm!=5)$V2  
  list2.perm <- subset(lhm.bim,V1_perm==5)$V2
  write.table(list1.perm,paste0(dir,"reml_XvsA_perm_circ_maf0.05/list_",paste0(i),"_1"),quote=F,row.names=F,col.names=F)
  write.table(list2.perm,paste0(dir,"reml_XvsA_perm_circ_maf0.05/list_",paste0(i),"_2"),quote=F,row.names=F,col.names=F)
}


#################################
## Male vs. female comparison
#################################

#Estimated Va, Autosome/X
fva <- read.table(paste0(dir,"reml_XvsA_maf0.05/XvsA.female.vars"),head=T)
mva <- read.table(paste0(dir,"reml_XvsA_maf0.05/XvsA.male.vars"),head=T)
#Adjust variances by a factor of two due to hemiclonal design
fva$Variance[fva$Component=="Var_K1"] <- 2*fva$Variance[fva$Component=="Var_K1"]
fva$SD[fva$Component=="Var_K1"] <- 2*fva$SD[fva$Component=="Var_K1"]
fva$Variance[fva$Component=="Var_K2"] <- 2*fva$Variance[fva$Component=="Var_K2"]
fva$SD[fva$Component=="Var_K2"] <- 2*fva$SD[fva$Component=="Var_K2"]
mva$Variance[mva$Component=="Var_K1"] <- 2*mva$Variance[mva$Component=="Var_K1"]
mva$SD[mva$Component=="Var_K1"] <- 2*mva$SD[mva$Component=="Var_K1"]

#Variances for each chromosome compartment
fva.A <- subset(fva,Component=="Var_K1")$Variance
fva.X <- subset(fva,Component=="Var_K2")$Variance
mva.A <- subset(mva,Component=="Var_K1")$Variance
mva.X <- subset(mva,Component=="Var_K2")$Variance

#Welch t-tests

#Autosome
t.test2(fva.A,mva.A,fva$SD[fva$Component=="Var_K1"]*sqrt(202),mva$SD[mva$Component=="Var_K1"]*sqrt(202),202,202,m0=0,equal.variance = F)
#p=0.019 

#X
t.test2(fva.X,mva.X,fva$SD[fva$Component=="Var_K2"]*sqrt(202),mva$SD[mva$Component=="Var_K2"]*sqrt(202),202,202,m0=0,equal.variance = F)
#p=0.676


#################################
## X vs Autosome comparison
#################################

#Circular permuted Va, for both sexes and both autosomes and X
fva.A.perm <- vector("numeric")
mva.A.perm <- vector("numeric")
fva.X.perm <- vector("numeric")
mva.X.perm <- vector("numeric")
for (i in 1:1000){
  fva.tmp <- read.table(paste0(dir,"reml_XvsA_perm_circ_maf0.05/reml_XvsA_perm/XvsA_female_",i,".vars"),head=T)
  fva.A.perm[i] <- 2*subset(fva.tmp,Component=="Var_K1")$Variance
  fva.X.perm[i] <- 2*subset(fva.tmp,Component=="Var_K2")$Variance
  mva.tmp <- read.table(paste0(dir,"reml_XvsA_perm_circ_maf0.05/reml_XvsA_perm/XvsA_male_",i,".vars"),head=T)
  mva.A.perm[i] <- 2*subset(mva.tmp,Component=="Var_K1")$Variance
  mva.X.perm[i] <- subset(mva.tmp,Component=="Var_K2")$Variance
}

#Rf
rf <- fva.X/(fva.A+fva.X)
rf.perm <- vector("numeric")
for (i in 1:1000){
  rf.perm[i] <- fva.X.perm[i]/(fva.A.perm[i]+fva.X.perm[i])
}
median(rf.perm)
#0.118
#p-value
sum(rf.perm<rf)/1000
#0.175

#Rm
rm <- mva.X/(mva.A+mva.X)
rm.perm <- vector("numeric")
for (i in 1:1000){
  rm.perm[i] <- mva.X.perm[i]/(mva.A.perm[i]+mva.X.perm[i])
}
median(rm.perm)
#0.048
#p-value
sum(rm.perm>rm)/1000
#0.174


#################################
## Plots
#################################

## Male vs female plot
mva$Sex <- "Male"
fva$Sex <- "Female"
bva <- rbind(fva,mva)
bva$Component <- c("Autosomes","X","E","Autosomes","X","E")
bva$Type <- "True"

pA <- ggplot(subset(bva,Component=="Autosomes" & Type=="True"),aes(y=Variance,x=Sex,col=Sex))+
  geom_errorbar(aes(ymin=Variance-SD,ymax=Variance+SD),width=0.3,position=position_dodge(0.2),size=0.8)+
  geom_point(position=position_dodge(0.2),size=2)+
  scale_color_manual(values=c("firebrick","dodgerblue4"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none",strip.text=element_text(size=15),strip.background=element_rect(fill="white"))+
  ylab(expression(paste("V"[A],"")))+
  geom_hline(yintercept=0)+
  ylim(c(-0.04,0.6))+
  xlab("")

pX <- ggplot(subset(bva,Component=="X" & Type=="True"),aes(y=Variance,x=Sex,col=Sex))+
  geom_errorbar(aes(ymin=Variance-SD,ymax=Variance+SD),width=0.3,position=position_dodge(0.2),size=0.8)+
  geom_point(position=position_dodge(0.2),size=2)+
  scale_color_manual(values=c("firebrick","dodgerblue4"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none",strip.text=element_text(size=15),strip.background=element_rect(fill="white"))+
  ylab(expression(paste("V"[A],"")))+
  geom_hline(yintercept=0)+
  ylim(c(-0.04,0.6))+
  xlab("")

p1 <- plot_grid(pA,pX,ncol = 2,labels = c("Autosomes","X chromosome"),scale=0.9,hjust = 0)
p1

## X vs autosome plot
bva.perm <- data.frame(matrix(data=NA,ncol=5,nrow=4000))
names(bva.perm) = names(bva)
bva.perm$Component <- c(rep("Autosomes",2000),rep("X",2000))
bva.perm$Variance <- c(fva.A.perm,mva.A.perm,fva.X.perm,mva.X.perm)
bva.perm$SD <- NA
bva.perm$Sex <- c(rep("Female",1000),rep("Male",1000),rep("Female",1000),rep("Male",1000))
bva.perm$Type <- "Permuted"
bva2 <- rbind(bva,bva.perm)

pF <- ggplot(subset(bva2,Component!="E" & Sex=="Female"),aes(y=Variance,x=Component))+
  geom_violin(data=subset(bva2,Type=="Permuted" & Sex=="Female"),width=0.7,size=1,col="grey75",fill="grey95")+
  geom_boxplot(data=subset(bva2,Type=="Permuted" & Sex=="Female"),width=0.1,alpha=0.5,size=1,col="grey75")+
  geom_point(data=subset(bva2,Component=="X" & Type=="True" & Sex=="Female"),shape=19,size=3,col="salmon")+
  geom_point(data=subset(bva2,Component=="Autosomes" & Type=="True" & Sex=="Female"),shape=19,size=3,col="forestgreen")+
  scale_size_manual(values=c(1,2))+
  scale_colour_manual(values=c("forestgreen",'salmon'))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.title=element_text(size=15),legend.text=element_text(size=15),strip.text=element_text(size=15),strip.background=element_rect(fill="white"))+
  xlab("")+
  ylab(expression(paste("V"[A],"")))+
  ylim(c(-0.04,0.6))+
  geom_hline(yintercept=0)

pM <- ggplot(subset(bva2,Component!="E" & Sex=="Male"),aes(y=Variance,x=Component))+
  geom_violin(data=subset(bva2,Type=="Permuted" & Sex=="Male"),width=0.7,size=1,col="grey75",fill="grey95")+
  geom_boxplot(data=subset(bva2,Type=="Permuted" & Sex=="Male"),width=0.1,alpha=0.5,size=1,col="grey75")+
  geom_point(data=subset(bva2,Component=="X" & Type=="True" & Sex=="Male"),shape=19,size=3,col="salmon")+
  geom_point(data=subset(bva2,Component=="Autosomes" & Type=="True" & Sex=="Male"),shape=19,size=3,col="forestgreen")+
  scale_size_manual(values=c(1,2))+
  scale_colour_manual(values=c("forestgreen",'salmon'))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.title=element_text(size=15),legend.text=element_text(size=15),strip.text=element_text(size=15),strip.background=element_rect(fill="white"))+
  xlab("")+
  ylab(expression(paste("V"[A],"")))+
  ylim(c(-0.04,0.6))+
  geom_hline(yintercept=0)

p2 <- plot_grid(pF,pM,ncol = 2,labels = c("Female","Male"),scale=0.9,hjust = 0)
p2

#Combine plots
alpha_minus_0.25 <- plot_grid(p1,p2,scale=0.9,ncol=2,labels = c("A","B"),hjust=0,rel_widths = c(0.5,0.5),label_size = 20)
alpha_minus_0.25

#################################
## Repeat with delta = -1
#################################

#Estimated Va, Autosome/X
fva <- read.table(paste0(dir,"reml_XvsA_maf0.05/XvsA.weight_minus_1.female.vars"),head=T)
mva <- read.table(paste0(dir,"reml_XvsA_maf0.05/XvsA.weight_minus_1.male.vars"),head=T)
#Adjust variances by a factor of two due to hemiclonal design
fva$Variance[fva$Component=="Var_K1"] <- 2*fva$Variance[fva$Component=="Var_K1"]
fva$SD[fva$Component=="Var_K1"] <- 2*fva$SD[fva$Component=="Var_K1"]
fva$Variance[fva$Component=="Var_K2"] <- 2*fva$Variance[fva$Component=="Var_K2"]
fva$SD[fva$Component=="Var_K2"] <- 2*fva$SD[fva$Component=="Var_K2"]
mva$Variance[mva$Component=="Var_K1"] <- 2*mva$Variance[mva$Component=="Var_K1"]
mva$SD[mva$Component=="Var_K1"] <- 2*mva$SD[mva$Component=="Var_K1"]

#Variances for each chromosome compartment
fva.A <- subset(fva,Component=="Var_K1")$Variance
fva.X <- subset(fva,Component=="Var_K2")$Variance
mva.A <- subset(mva,Component=="Var_K1")$Variance
mva.X <- subset(mva,Component=="Var_K2")$Variance

#Welch t-tests

#Autosome
t.test2(fva.A,mva.A,fva$SD[fva$Component=="Var_K1"]*sqrt(202),mva$SD[mva$Component=="Var_K1"]*sqrt(202),202,202,m0=0,equal.variance = F)
#p=0.029 

#X
t.test2(fva.X,mva.X,fva$SD[fva$Component=="Var_K2"]*sqrt(202),mva$SD[mva$Component=="Var_K2"]*sqrt(202),202,202,m0=0,equal.variance = F)
#p=0.772


#Circular permuted Va, for both sexes and both autosomes and X
fva.A.perm <- vector("numeric")
mva.A.perm <- vector("numeric")
fva.X.perm <- vector("numeric")
mva.X.perm <- vector("numeric")
for (i in 1:1000){
  fva.tmp <- read.table(paste0(dir,"reml_XvsA_perm_circ_maf0.05/reml_XvsA_perm_weight_minus_1/XvsA_female_",i,".vars"),head=T)
  fva.A.perm[i] <- 2*subset(fva.tmp,Component=="Var_K1")$Variance
  fva.X.perm[i] <- 2*subset(fva.tmp,Component=="Var_K2")$Variance
  mva.tmp <- read.table(paste0(dir,"reml_XvsA_perm_circ_maf0.05/reml_XvsA_perm_weight_minus_1/XvsA_male_",i,".vars"),head=T)
  mva.A.perm[i] <- 2*subset(mva.tmp,Component=="Var_K1")$Variance
  mva.X.perm[i] <- subset(mva.tmp,Component=="Var_K2")$Variance
}


#Rf
rf <- fva.X/(fva.A+fva.X)
rf.perm <- vector("numeric")
for (i in 1:1000){
  rf.perm[i] <- fva.X.perm[i]/(fva.A.perm[i]+fva.X.perm[i]) 
}
median(rf.perm)
#0.118
#p-value
sum(rf.perm<rf)/1000
#0.279

#Rm
rm <-  mva.X/(mva.A+mva.X) 
rm.perm <- vector("numeric")
for (i in 1:1000){
  rm.perm[i] <- mva.X.perm[i]/(mva.A.perm[i]+mva.X.perm[i]) 
}
median(rm.perm)
#0.048
#p-value
sum(rm.perm>rm)/1000
#0.173


## Male vs female plot
mva$Sex <- "Male"
fva$Sex <- "Female"
bva <- rbind(fva,mva)
bva$Component <- c("Autosomes","X","E","Autosomes","X","E")
bva$Type <- "True"

pA <- ggplot(subset(bva,Component=="Autosomes" & Type=="True"),aes(y=Variance,x=Sex,col=Sex))+
  geom_errorbar(aes(ymin=Variance-SD,ymax=Variance+SD),width=0.3,position=position_dodge(0.2),size=0.8)+
  geom_point(position=position_dodge(0.2),size=2)+
  scale_color_manual(values=c("firebrick","dodgerblue4"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none",strip.text=element_text(size=15),strip.background=element_rect(fill="white"))+
  ylab(expression(paste("V"[A],"")))+
  geom_hline(yintercept=0)+
  ylim(c(-0.04,0.6))+
  xlab("")

pX <- ggplot(subset(bva,Component=="X" & Type=="True"),aes(y=Variance,x=Sex,col=Sex))+
  geom_errorbar(aes(ymin=Variance-SD,ymax=Variance+SD),width=0.3,position=position_dodge(0.2),size=0.8)+
  geom_point(position=position_dodge(0.2),size=2)+
  scale_color_manual(values=c("firebrick","dodgerblue4"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none",strip.text=element_text(size=15),strip.background=element_rect(fill="white"))+
  ylab(expression(paste("V"[A],"")))+
  geom_hline(yintercept=0)+
  ylim(c(-0.04,0.6))+
  xlab("")

p1 <- plot_grid(pA,pX,ncol = 2,labels = c("Autosomes","X chromosome"),scale=0.9,hjust = 0)
p1

## X vs autosome plot
bva.perm <- data.frame(matrix(data=NA,ncol=5,nrow=4000))
names(bva.perm) = names(bva)
bva.perm$Component <- c(rep("Autosomes",2000),rep("X",2000))
bva.perm$Variance <- c(fva.A.perm,mva.A.perm,fva.X.perm,mva.X.perm)
bva.perm$SD <- NA
bva.perm$Sex <- c(rep("Female",1000),rep("Male",1000),rep("Female",1000),rep("Male",1000))
bva.perm$Type <- "Permuted"
bva2 <- rbind(bva,bva.perm)

pF <- ggplot(subset(bva2,Component!="E" & Sex=="Female"),aes(y=Variance,x=Component))+
  geom_violin(data=subset(bva2,Type=="Permuted" & Sex=="Female"),width=0.7,size=1,col="grey75",fill="grey95")+
  geom_boxplot(data=subset(bva2,Type=="Permuted" & Sex=="Female"),width=0.1,alpha=0.5,size=1,col="grey75")+
  geom_point(data=subset(bva2,Component=="X" & Type=="True" & Sex=="Female"),shape=19,size=3,col="salmon")+
  geom_point(data=subset(bva2,Component=="Autosomes" & Type=="True" & Sex=="Female"),shape=19,size=3,col="forestgreen")+
  scale_size_manual(values=c(1,2))+
  scale_colour_manual(values=c("forestgreen",'salmon'))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.title=element_text(size=15),legend.text=element_text(size=15),strip.text=element_text(size=15),strip.background=element_rect(fill="white"))+
  xlab("")+
  ylab(expression(paste("V"[A],"")))+
  ylim(c(-0.04,0.6))+
  geom_hline(yintercept=0)

pM <- ggplot(subset(bva2,Component!="E" & Sex=="Male"),aes(y=Variance,x=Component))+
  geom_violin(data=subset(bva2,Type=="Permuted" & Sex=="Male"),width=0.7,size=1,col="grey75",fill="grey95")+
  geom_boxplot(data=subset(bva2,Type=="Permuted" & Sex=="Male"),width=0.1,alpha=0.5,size=1,col="grey75")+
  geom_point(data=subset(bva2,Component=="X" & Type=="True" & Sex=="Male"),shape=19,size=3,col="salmon")+
  geom_point(data=subset(bva2,Component=="Autosomes" & Type=="True" & Sex=="Male"),shape=19,size=3,col="forestgreen")+
  scale_size_manual(values=c(1,2))+
  scale_colour_manual(values=c("forestgreen",'salmon'))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.title=element_text(size=15),legend.text=element_text(size=15),strip.text=element_text(size=15),strip.background=element_rect(fill="white"))+
  xlab("")+
  ylab(expression(paste("V"[A],"")))+
  ylim(c(-0.04,0.6))+
  geom_hline(yintercept=0)

p2 <- plot_grid(pF,pM,ncol = 2,labels = c("Female","Male"),scale=0.9,hjust = 0)
p2

#Combine plots
alpha_minus_1 <- plot_grid(p1,p2,scale=0.9,ncol=2,labels = c("",""),hjust=0,rel_widths = c(0.5,0.5),label_size = 20)



#################################
## Repeat with delta = 0
#################################


#Estimated Va, Autosome/X
fva <- read.table(paste0(dir,"reml_XvsA_maf0.05/XvsA.weight_0.female.vars"),head=T)
mva <- read.table(paste0(dir,"reml_XvsA_maf0.05/XvsA.weight_0.male.vars"),head=T)
#Adjust variances by a factor of two due to hemiclonal design
fva$Variance[fva$Component=="Var_K1"] <- 2*fva$Variance[fva$Component=="Var_K1"]
fva$SD[fva$Component=="Var_K1"] <- 2*fva$SD[fva$Component=="Var_K1"]
fva$Variance[fva$Component=="Var_K2"] <- 2*fva$Variance[fva$Component=="Var_K2"]
fva$SD[fva$Component=="Var_K2"] <- 2*fva$SD[fva$Component=="Var_K2"]
mva$Variance[mva$Component=="Var_K1"] <- 2*mva$Variance[mva$Component=="Var_K1"]
mva$SD[mva$Component=="Var_K1"] <- 2*mva$SD[mva$Component=="Var_K1"]

#Variances for each chromosome compartment
fva.A <- subset(fva,Component=="Var_K1")$Variance
fva.X <- subset(fva,Component=="Var_K2")$Variance
mva.A <- subset(mva,Component=="Var_K1")$Variance
mva.X <- subset(mva,Component=="Var_K2")$Variance

#Welch t-tests

#Autosome
t.test2(fva.A,mva.A,fva$SD[fva$Component=="Var_K1"]*sqrt(202),mva$SD[mva$Component=="Var_K1"]*sqrt(202),202,202,m0=0,equal.variance = F)
#p=0.018

#X
t.test2(fva.X,mva.X,fva$SD[fva$Component=="Var_K2"]*sqrt(202),mva$SD[mva$Component=="Var_K2"]*sqrt(202),202,202,m0=0,equal.variance = F)
#p=0.658


#Circular permuted Va, for both sexes and both autosomes and X
fva.A.perm <- vector("numeric")
mva.A.perm <- vector("numeric")
fva.X.perm <- vector("numeric")
mva.X.perm <- vector("numeric")
for (i in 1:1000){
  fva.tmp <- read.table(paste0(dir,"reml_XvsA_perm_circ_maf0.05/reml_XvsA_perm_weight_0/XvsA_female_",i,".vars"),head=T)
  fva.A.perm[i] <- 2*subset(fva.tmp,Component=="Var_K1")$Variance
  fva.X.perm[i] <- 2*subset(fva.tmp,Component=="Var_K2")$Variance
  mva.tmp <- read.table(paste0(dir,"reml_XvsA_perm_circ_maf0.05/reml_XvsA_perm_weight_0/XvsA_male_",i,".vars"),head=T)
  mva.A.perm[i] <- 2*subset(mva.tmp,Component=="Var_K1")$Variance
  mva.X.perm[i] <- subset(mva.tmp,Component=="Var_K2")$Variance
}


#Rf
rf <- fva.X/(fva.A+fva.X) 
rf.perm <- vector("numeric")
for (i in 1:1000){
  rf.perm[i] <- fva.X.perm[i]/(fva.A.perm[i]+fva.X.perm[i]) 
}
median(rf.perm)
#0.116
#p-value
sum(rf.perm<rf)/1000
#0.175

#Rm
rm <- mva.X/(mva.A+mva.X) 
rm.perm <- vector("numeric")
for (i in 1:1000){
  rm.perm[i] <- mva.X.perm[i]/(mva.A.perm[i]+mva.X.perm[i]) 
}
median(rm.perm)
#0.044
#p-value
sum(rm.perm>rm)/1000
#0.174


## Male vs female plot
mva$Sex <- "Male"
fva$Sex <- "Female"
bva <- rbind(fva,mva)
bva$Component <- c("Autosomes","X","E","Autosomes","X","E")
bva$Type <- "True"

pA <- ggplot(subset(bva,Component=="Autosomes" & Type=="True"),aes(y=Variance,x=Sex,col=Sex))+
  geom_errorbar(aes(ymin=Variance-SD,ymax=Variance+SD),width=0.3,position=position_dodge(0.2),size=0.8)+
  geom_point(position=position_dodge(0.2),size=2)+
  scale_color_manual(values=c("firebrick","dodgerblue4"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none",strip.text=element_text(size=15),strip.background=element_rect(fill="white"))+
  ylab(expression(paste("V"[A],"")))+
  geom_hline(yintercept=0)+
  ylim(c(-0.04,0.6))+
  xlab("")

pX <- ggplot(subset(bva,Component=="X" & Type=="True"),aes(y=Variance,x=Sex,col=Sex))+
  geom_errorbar(aes(ymin=Variance-SD,ymax=Variance+SD),width=0.3,position=position_dodge(0.2),size=0.8)+
  geom_point(position=position_dodge(0.2),size=2)+
  scale_color_manual(values=c("firebrick","dodgerblue4"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none",strip.text=element_text(size=15),strip.background=element_rect(fill="white"))+
  ylab(expression(paste("V"[A],"")))+
  geom_hline(yintercept=0)+
  ylim(c(-0.04,0.6))+
  xlab("")

p1 <- plot_grid(pA,pX,ncol = 2,labels = c("Autosomes","X chromosome"),scale=0.9,hjust = 0)
p1

## X vs autosome plot
bva.perm <- data.frame(matrix(data=NA,ncol=5,nrow=4000))
names(bva.perm) = names(bva)
bva.perm$Component <- c(rep("Autosomes",2000),rep("X",2000))
bva.perm$Variance <- c(fva.A.perm,mva.A.perm,fva.X.perm,mva.X.perm)
bva.perm$SD <- NA
bva.perm$Sex <- c(rep("Female",1000),rep("Male",1000),rep("Female",1000),rep("Male",1000))
bva.perm$Type <- "Permuted"
bva2 <- rbind(bva,bva.perm)

pF <- ggplot(subset(bva2,Component!="E" & Sex=="Female"),aes(y=Variance,x=Component))+
  geom_violin(data=subset(bva2,Type=="Permuted" & Sex=="Female"),width=0.7,size=1,col="grey75",fill="grey95")+
  geom_boxplot(data=subset(bva2,Type=="Permuted" & Sex=="Female"),width=0.1,alpha=0.5,size=1,col="grey75")+
  geom_point(data=subset(bva2,Component=="X" & Type=="True" & Sex=="Female"),shape=19,size=3,col="salmon")+
  geom_point(data=subset(bva2,Component=="Autosomes" & Type=="True" & Sex=="Female"),shape=19,size=3,col="forestgreen")+
  scale_size_manual(values=c(1,2))+
  scale_colour_manual(values=c("forestgreen",'salmon'))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.title=element_text(size=15),legend.text=element_text(size=15),strip.text=element_text(size=15),strip.background=element_rect(fill="white"))+
  xlab("")+
  ylab(expression(paste("V"[A],"")))+
  ylim(c(-0.04,0.6))+
  geom_hline(yintercept=0)

pM <- ggplot(subset(bva2,Component!="E" & Sex=="Male"),aes(y=Variance,x=Component))+
  geom_violin(data=subset(bva2,Type=="Permuted" & Sex=="Male"),width=0.7,size=1,col="grey75",fill="grey95")+
  geom_boxplot(data=subset(bva2,Type=="Permuted" & Sex=="Male"),width=0.1,alpha=0.5,size=1,col="grey75")+
  geom_point(data=subset(bva2,Component=="X" & Type=="True" & Sex=="Male"),shape=19,size=3,col="salmon")+
  geom_point(data=subset(bva2,Component=="Autosomes" & Type=="True" & Sex=="Male"),shape=19,size=3,col="forestgreen")+
  scale_size_manual(values=c(1,2))+
  scale_colour_manual(values=c("forestgreen",'salmon'))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.title=element_text(size=15),legend.text=element_text(size=15),strip.text=element_text(size=15),strip.background=element_rect(fill="white"))+
  xlab("")+
  ylab(expression(paste("V"[A],"")))+
  ylim(c(-0.04,0.6))+
  geom_hline(yintercept=0)

p2 <- plot_grid(pF,pM,ncol = 2,labels = c("Female","Male"),scale=0.9,hjust = 0)
p2

alpha_minus_0 <- plot_grid(p1,p2,scale=0.9,ncol=2,labels = c("",""),hjust=0,rel_widths = c(0.5,0.5),label_size = 20)

#Combine plots^2
plot_grid(alpha_minus_1,alpha_minus_0,scale=0.9,ncol=1,hjust=0,labels = c("delta = -1","delta = 0"),label_size = 15)

