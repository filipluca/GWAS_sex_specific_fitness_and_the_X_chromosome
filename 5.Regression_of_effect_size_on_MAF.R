library(CMplot)
library(qqman)
library(Hmisc)
library(reshape2)
library(ggplot2)
library(data.table)
library(Rmisc)
library(cowplot)

dir <- "~/Dropbox/dropbox_work/data/GWAS_m_f_fitness/"

################################################
##Effect size ~ MAF####
##Male vs. female comparison ####
################################################

## Import unpermuted GWAS data
fassoc <- read.table(paste0(dir,"female_maf0.05.assoc"),head=T)
massoc <- read.table(paste0(dir,"male_maf0.05.assoc"),head=T)
#Absolute effect size
fassoc$Abs_Effect <- abs(fassoc$Effect)
massoc$Abs_Effect <- abs(massoc$Effect)
#Chromosome coding
fassoc$Chromosome <- as.factor(fassoc$Chromosome)
levels(fassoc$Chromosome) <- c("2L","2R","3L","3R","X")
massoc$Chromosome <- as.factor(massoc$Chromosome)
levels(massoc$Chromosome) <- c("2L","2R","3L","3R","X")
#Title of predictor column from 'Predictor' -> 'SNP'
names(fassoc)[2] <- "SNP"
names(massoc)[2] <- "SNP"
#LD pruned predictors
fassoc.pruned <- read.table(paste0(dir,"female_maf0.05_pruned_r0.4_10kb.clumped"),head=T)
fassoc$LD_indep <- ifelse(fassoc$SNP %in% fassoc.pruned$SNP,1,0)
massoc$LD_indep <- ifelse(massoc$SNP %in% fassoc.pruned$SNP,1,0)

## Permuted phenotypes
female.file.names <- list.files(paste0(dir,"1000_permuted_phenos/"),pattern = "^female_*")
male.file.names <- list.files(paste0(dir,"1000_permuted_phenos/"),pattern = "^male_*")

#GLM coefficients for...
#Sex * MAF (autosomes)
Gamma_GLM_effects_A <- matrix(data=NA,ncol = 3,nrow=1000)
#Sex * MAF (X)
Gamma_GLM_effects_X <- matrix(data=NA,ncol = 3,nrow=1000)
#Chromosome * MAF (males)
Gamma_GLM_effects_m <- matrix(data=NA,ncol = 3,nrow=1000)
#Chromosome * MAF (females)
Gamma_GLM_effects_f <- matrix(data=NA,ncol = 3,nrow=1000)
set.seed(123)
for (i in 1:1000){
  #Import Female
  fassoc.perm <- fread(paste0(dir,"1000_permuted_phenos/",female.file.names[i]),head=T)
  fassoc.perm$Abs_Effect <- abs(fassoc.perm$BETA)
  fassoc.perm$Chromosome <- fassoc$Chromosome
  fassoc.perm$MAF <- fassoc$MAF
  fassoc.perm$LD_indep <- fassoc$LD_indep
  fassoc.perm$Sex <- "Female"
  #Import Male
  massoc.perm <- fread(paste0(dir,"1000_permuted_phenos/",male.file.names[i]),head=T)
  massoc.perm$Abs_Effect <- abs(massoc.perm$BETA)
  massoc.perm$Chromosome <- massoc$Chromosome
  massoc.perm$MAF <- massoc$MAF
  massoc.perm$LD_indep <- massoc$LD_indep
  massoc.perm$Sex <- "Male"
  
  #Gamma GLM effects
  bassoc.perm <- rbind(fassoc.perm,massoc.perm)
  Gamma_GLM_effects_A[i,1:3] <- as.numeric(coef(glm(data=subset(bassoc.perm,Abs_Effect>0 & LD_indep==1 & Chromosome!="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[2:4])
  Gamma_GLM_effects_X[i,1:3] <- as.numeric(coef(glm(data=subset(bassoc.perm,Abs_Effect>0 & LD_indep==1 & Chromosome=="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[2:4])

  Gamma_GLM_effects_m[i,1:3] <- as.numeric(coef(glm(data=subset(massoc.perm,Abs_Effect>0 & LD_indep==1 ),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[2:4])
  Gamma_GLM_effects_f[i,1:3] <- as.numeric(coef(glm(data=subset(fassoc.perm,Abs_Effect>0 & LD_indep==1 ),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[2:4])
  print(i)
  }

##P-values for each model term (autosomes)
fassoc$Sex <- "Female"
massoc$Sex <- "Male"
bassoc <- rbind(fassoc,massoc)

#Is Gamma glm appropriate?
model <- glm(data=subset(bassoc,Abs_Effect>0 & LD_indep==1 & Chromosome!="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log"))
model <- glm(data=subset(bassoc,Abs_Effect>0 & LD_indep==1 & Chromosome=="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log"))
model <- glm(data=subset(massoc,Abs_Effect>0 & LD_indep==1 ),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log"))
model <- glm(data=subset(fassoc,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log"))
ypred <-  predict(model)
res <-  residuals(model, type = 'deviance')
plot(ypred,res)
hist(res)

#Sex effect
sum(coef(glm(data=subset(bassoc,Abs_Effect>0 & LD_indep==1 & Chromosome!="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[2]>Gamma_GLM_effects_A[,1])
#p<0.001 (p<0.001 using non-squared abs effect)
#MAF effect
sum(coef(glm(data=subset(bassoc,Abs_Effect>0 & LD_indep==1 & Chromosome!="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[3]>Gamma_GLM_effects_A[,2])
#p=0.830 (p=0.732 using non-squared abs effect)
#Sex*MAF effect
sum(coef(glm(data=subset(bassoc,Abs_Effect>0 & LD_indep==1 & Chromosome!="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[4]<Gamma_GLM_effects_A[,3])
#p=0.287 (p=0.205 using non-squared abs effect)

##P-values for each model term (X)
#Sex effect
sum(coef(glm(data=subset(bassoc,Abs_Effect>0 & LD_indep==1 & Chromosome=="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[2]>Gamma_GLM_effects_X[,1])
#p=0.001 (p<0.001 using non-squared abs effect)
#MAF effect
sum(coef(glm(data=subset(bassoc,Abs_Effect>0 & LD_indep==1 & Chromosome=="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[3]>Gamma_GLM_effects_X[,2])
#p=0.162 (p=0.104 using non-squared abs effect)
#Sex*MAF effect
sum(coef(glm(data=subset(bassoc,Abs_Effect>0 & LD_indep==1 & Chromosome=="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[4]<Gamma_GLM_effects_X[,3])
#p=0.429 (p=0.207 using non-squared abs effect)

## ES~MAF plot (autosomes)
p1 <- ggplot(subset(bassoc,Abs_Effect>0 & LD_indep==1& Chromosome!="X"),aes(y=Abs_Effect,x=MAF,col=Sex))+
  geom_point(alpha=0.4)+
  scale_color_manual(values=c("firebrick","dodgerblue4"))+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")
p2 <- ggplot(subset(bassoc,Abs_Effect>0 & LD_indep==1& Chromosome!="X"),aes(y=Abs_Effect,x=MAF,col=Sex))+
  scale_color_manual(values=c("firebrick","dodgerblue4"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")+
  geom_smooth(alpha=1,method = "glm", method.args = list(family = "Gamma"))
mvsf.A <- ggdraw()+draw_plot(p1)+draw_plot(p2,x = 0.57, y = 0.57, width = 2/5, height =  2/5)

## ES~MAF plot (X chromosome)
p1 <- ggplot(subset(bassoc,Abs_Effect>0 & LD_indep==1& Chromosome=="X"),aes(y=Abs_Effect,x=MAF,col=Sex))+
  geom_point(alpha=0.4)+
  scale_color_manual(values=c("firebrick","dodgerblue4"))+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")
p2 <- ggplot(subset(bassoc,Abs_Effect>0 & LD_indep==1& Chromosome=="X"),aes(y=Abs_Effect,x=MAF,col=Sex))+
  scale_color_manual(values=c("firebrick","dodgerblue4"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")+
  geom_smooth(alpha=1,method = "glm", method.args = list(family = "Gamma"))
mvsf.X <- ggdraw()+draw_plot(p1)+draw_plot(p2,x = 0.57, y = 0.57, width = 2/5, height =  2/5)

## SupFig (scatter plots of observed and permuted model terms)

#Autosomes
Gamma_GLM_effects_A_df <- as.data.frame(Gamma_GLM_effects_A)
names(Gamma_GLM_effects_A_df) <- c("Sex","MAF","Sex_by_MAF")
Gamma_GLM_effects_A_df_melted <- melt(Gamma_GLM_effects_A_df)
Gamma_GLM_effects_A_df_melted$type <- "Permuted"
Gamma_GLM_effects_A_df_melted[3001:3003,1] <- c("Sex","MAF","Sex_by_MAF")
Gamma_GLM_effects_A_df_melted[3001:3003,2] <- coef(glm(data=subset(bassoc,Abs_Effect>0 & LD_indep==1 & Chromosome!="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[2:4]
Gamma_GLM_effects_A_df_melted[3001:3003,3] <- rep("True",3)
levels(Gamma_GLM_effects_A_df_melted$variable) <- c("Sex (male)","MAF","Sex(male)-by-MAF")

#X 
Gamma_GLM_effects_X_df <- as.data.frame(Gamma_GLM_effects_X)
names(Gamma_GLM_effects_X_df) <- c("Sex","MAF","Sex_by_MAF")
Gamma_GLM_effects_X_df_melted <- melt(Gamma_GLM_effects_X_df)
Gamma_GLM_effects_X_df_melted$type <- "Permuted"
Gamma_GLM_effects_X_df_melted[3001:3003,1] <- c("Sex","MAF","Sex_by_MAF")
Gamma_GLM_effects_X_df_melted[3001:3003,2] <- coef(glm(data=subset(bassoc,Abs_Effect>0 & LD_indep==1 & Chromosome=="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[2:4]
Gamma_GLM_effects_X_df_melted[3001:3003,3] <- rep("True",3)
levels(Gamma_GLM_effects_X_df_melted$variable) <- c("Sex (male)","MAF","Sex(male)-by-MAF")

#Combine X and autosome 
Gamma_GLM_effects_A_df_melted$Chromosome <- "Autosomes"
Gamma_GLM_effects_X_df_melted$Chromosome <- "X"
Gamma_GLM_effects_df_melted <- rbind(Gamma_GLM_effects_A_df_melted,Gamma_GLM_effects_X_df_melted)
pApX <- ggplot(Gamma_GLM_effects_df_melted,aes(x=Chromosome,y=value,size=type))+
  geom_jitter(data=subset(Gamma_GLM_effects_df_melted,type=="Permuted" & variable=="Sex(male)-by-MAF"),width=0.1,alpha=0.5,size=1.5,col="grey")+
  geom_point(data=subset(Gamma_GLM_effects_df_melted,type=="True"& variable=="Sex(male)-by-MAF"),shape=19,size=3,col="black")+
  #geom_point(data=subset(Gamma_GLM_effects_X_df_melted,type=="True" & variable=="Sex (male)"),shape=19,size=3,col="dodgerblue4")+
  scale_size_manual(values=c(1,2))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none")+
  xlab("")+
  ylab("Regression coefficient \n(sex-by-MAF interaction)")


################################################
##Effect size ~ MAF####
##X vs. autosome comparison ####
################################################

## Permuted genotypic data 
Gamma_GLM_effects_circ_f <- matrix(data=NA,ncol = 3,nrow=1000)
Gamma_GLM_effects_circ_m <- matrix(data=NA,ncol = 3,nrow=1000)
start.point <- vector("numeric")
set.seed(123)
for (i in 1:1000){
  #Initiate X chromosome at a random starting point
  chrom.circular <- rep(fassoc$Chromosome,2)
  start.point[i] <- sample(1:765764,1)
  fassoc$Chromosome_perm <- chrom.circular[start.point[i]:(start.point[i]+765763)]
  massoc$Chromosome_perm <- chrom.circular[start.point[i]:(start.point[i]+765763)]
  
  #Gamma GLM effects
  Gamma_GLM_effects_circ_f[i,1:3] <- as.numeric(coef(glm(data=subset(fassoc,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome_perm=="X")*MAF,family=Gamma(link="log")))[2:4])
  Gamma_GLM_effects_circ_m[i,1:3] <- as.numeric(coef(glm(data=subset(massoc,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome_perm=="X")*MAF,family=Gamma(link="log")))[2:4])
  print(i)
}

## P-values for each model term, autosomes vs X (females)
#Chromosome effect
Chr.effect <- as.numeric(coef(glm(data=subset(fassoc,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[2])
sum(Chr.effect>Gamma_GLM_effects_circ_f[,1])/1000
#p=0.381 (p=0.451, using non-squared abs effect)
#MAF effect
MAF.effect <- as.numeric(coef(glm(data=subset(fassoc,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[3])
sum(MAF.effect>Gamma_GLM_effects_circ_f[,2])/1000
#p=0.936 (p=0.935, using non-squared abs effect)
#Interaction effect
X.int.effect <- as.numeric(coef(glm(data=subset(fassoc,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[4])
sum(X.int.effect>Gamma_GLM_effects_circ_f[,3])/1000
#p=0.033 (p=0.048, using non-squared abs effect)

## P-values for each model term, autosomes vs X (males)
#Chromosome effect
Chr.effect <- as.numeric(coef(glm(data=subset(massoc,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[2])
sum(Chr.effect>Gamma_GLM_effects_circ_m[,1])/1000
#p=0.854 (p=0.832, using non-squared abs effect)
#MAF effect
MAF.effect <- as.numeric(coef(glm(data=subset(massoc,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[3])
sum(MAF.effect>Gamma_GLM_effects_circ_m[,2])/1000
#p=0.933 (p=0.925, using non-squared abs effects)
#Interaction effect
X.int.effect <- as.numeric(coef(glm(data=subset(massoc,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[4])
sum(X.int.effect>Gamma_GLM_effects_circ_m[,3])/1000
#p=0.095 (p=0.140, using non-squared abs effects)


## ES~MAF plot (females)
fassoc$Chromosome_binary <- ifelse(fassoc$Chromosome=='X',"X","Autosomes")
p1 <- ggplot(subset(fassoc,Abs_Effect>0 & LD_indep==1),aes(y=Abs_Effect,x=MAF,col=Chromosome_binary))+
  geom_point(alpha=0.4)+
  scale_color_manual(values=c("forestgreen","salmon"))+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")
p2 <- ggplot(subset(fassoc,Abs_Effect>0 & LD_indep==1),aes(y=Abs_Effect,x=MAF,col=Chromosome_binary))+
  scale_color_manual(values=c("forestgreen","salmon"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")+
  geom_smooth(alpha=1,method = "glm", method.args = list(family = "Gamma"))
XvsA.f <- ggdraw()+draw_plot(p1)+draw_plot(p2,x = 0.57, y = 0.57, width = 2/5, height =  2/5)

## ES~MAF plot (males)
massoc$Chromosome_binary <- ifelse(massoc$Chromosome=='X',"X","Autosomes")
p1 <- ggplot(subset(massoc,Abs_Effect>0 & LD_indep==1),aes(y=Abs_Effect,x=MAF,col=Chromosome_binary))+
  geom_point(alpha=0.4)+
  scale_color_manual(values=c("forestgreen","salmon"))+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")
p2 <- ggplot(subset(massoc,Abs_Effect>0 & LD_indep==1),aes(y=Abs_Effect,x=MAF,col=Chromosome_binary))+
  scale_color_manual(values=c("forestgreen","salmon"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")+
  geom_smooth(alpha=1,method = "glm", method.args = list(family = "Gamma"))
XvsA.m <- ggdraw()+draw_plot(p1)+draw_plot(p2,x = 0.57, y = 0.57, width = 2/5, height =  2/5)

## SupFig (scatter plot of model coefficients)
#Females
Gamma_GLM_effects_circ_f <- as.data.frame(Gamma_GLM_effects_circ_f)
names(Gamma_GLM_effects_circ_f)[1:3] <- c("Chrom.(X)","MAF","Chrom.(X)-by-MAF")
Gamma_GLM_effects_circ_f$Mode <- "Genotype"
Gamma_GLM_effects_f_melted1 <- melt(Gamma_GLM_effects_circ_f)
Gamma_GLM_effects_f <- as.data.frame(Gamma_GLM_effects_f)
names(Gamma_GLM_effects_f)[1:3] <- c("Chrom.(X)","MAF","Chrom.(X)-by-MAF")
Gamma_GLM_effects_f$Mode <- "Phenotype"
Gamma_GLM_effects_f_melted2 <- melt(Gamma_GLM_effects_f)
Gamma_GLM_effects_f_melted <- rbind(Gamma_GLM_effects_f_melted1,Gamma_GLM_effects_f_melted2)
Gamma_GLM_effects_f_melted$type <- "Permuted"
Gamma_GLM_effects_f_melted[6001:6003,2] <- c("Chrom.(X)","MAF","Chrom.(X)-by-MAF")
Gamma_GLM_effects_f_melted[6001:6003,3] <- coef(glm(data=subset(fassoc,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[2:4]
Gamma_GLM_effects_f_melted[6001:6003,4] <- rep("True",3)
Gamma_GLM_effects_f_melted[6001:6003,1] <- NA
Gamma_GLM_effects_f_melted$value <- as.numeric(as.character(Gamma_GLM_effects_f_melted$value))

#Males
Gamma_GLM_effects_circ_m <- as.data.frame(Gamma_GLM_effects_circ_m)
names(Gamma_GLM_effects_circ_m)[1:3] <- c("Chrom.(X)","MAF","Chrom.(X)-by-MAF")
Gamma_GLM_effects_circ_m$Mode <- "Genotype"
Gamma_GLM_effects_m_melted1 <- melt(Gamma_GLM_effects_circ_m)
Gamma_GLM_effects_m <- as.data.frame(Gamma_GLM_effects_m)
names(Gamma_GLM_effects_m)[1:3] <- c("Chrom.(X)","MAF","Chrom.(X)-by-MAF")
Gamma_GLM_effects_m$Mode <- "Phenotype"
Gamma_GLM_effects_m_melted2 <- melt(Gamma_GLM_effects_m)
Gamma_GLM_effects_m_melted <- rbind(Gamma_GLM_effects_m_melted1,Gamma_GLM_effects_m_melted2)
Gamma_GLM_effects_m_melted$type <- "Permuted"
Gamma_GLM_effects_m_melted[6001:6003,2] <- c("Chrom.(X)","MAF","Chrom.(X)-by-MAF")
Gamma_GLM_effects_m_melted[6001:6003,3] <- coef(glm(data=subset(massoc,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[2:4]
Gamma_GLM_effects_m_melted[6001:6003,4] <- rep("True",3)
Gamma_GLM_effects_m_melted[6001:6003,1] <- NA
Gamma_GLM_effects_m_melted$value <- as.numeric(as.character(Gamma_GLM_effects_m_melted$value))

#Combine females and males
Gamma_GLM_effects_f_melted$Sex <- "Female"
Gamma_GLM_effects_m_melted$Sex <- "Male"
Gamma_GLM_effects_fm_melted <- rbind(Gamma_GLM_effects_f_melted,Gamma_GLM_effects_m_melted)
pFpM <- ggplot(Gamma_GLM_effects_fm_melted,aes(x=Sex,y=value,col=type,size=type))+
  #geom_jitter(data=subset(Gamma_GLM_effects_f_melted,Mode=="Phenotype"),width=0.1,alpha=0.5,size=1.5,col="grey")+
  geom_jitter(data=subset(Gamma_GLM_effects_fm_melted,Mode=="Genotype" & variable=="Chrom.(X)-by-MAF"),width=0.1,alpha=0.5,size=1.5,col="grey")+
  geom_point(data=subset(Gamma_GLM_effects_fm_melted,type=="True" & variable=="Chrom.(X)-by-MAF"),shape=19,size=3,col="black")+
  #geom_point(data=subset(Gamma_GLM_effects_fm_melted,type=="True" & variable=="Chrom.(X)"),shape=19,size=3,col="salmon")+
  scale_size_manual(values=c(1,2))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none")+
  xlab("")+
  ylab("Regression coefficient \n(Compartment-by-MAF interaction)")


################################################
## Delta = -1####
## Male vs. female####
################################################

## Import GWAS data
fassoc1 <- read.table(paste0(dir,"female_maf0.05_weight_minus_1.assoc"),head=T)
massoc1 <- read.table(paste0(dir,"male_maf0.05_weight_minus_1.assoc"),head=T)
#Absolute effect size
fassoc1$Abs_Effect <- abs(fassoc1$Effect)
massoc1$Abs_Effect <- abs(massoc1$Effect)
#Chromosome coding
fassoc1$Chromosome <- as.factor(fassoc1$Chromosome)
levels(fassoc1$Chromosome) <- c("2L","2R","3L","3R","X")
massoc1$Chromosome <- as.factor(massoc1$Chromosome)
levels(massoc1$Chromosome) <- c("2L","2R","3L","3R","X")
#Title of predictor column from 'Predictor' -> 'SNP' (for LD pruning)
names(fassoc1)[2] <- "SNP"
names(massoc1)[2] <- "SNP"
#LD pruned predictors
fassoc.pruned <- read.table(paste0(dir,"female_maf0.05_pruned_r0.4_10kb.clumped"),head=T)
fassoc1$LD_indep <- ifelse(fassoc1$SNP %in% fassoc.pruned$SNP,1,0)
massoc1$LD_indep <- ifelse(massoc1$SNP %in% fassoc.pruned$SNP,1,0)

## Permuted phenos
female.file.names <- list.files(paste0(dir,"1000_permuted_phenos_weight_minus_1/"),pattern = "^female_*")
male.file.names <- list.files(paste0(dir,"1000_permuted_phenos_weight_minus_1/"),pattern = "^male_*")

#GLM coefficients for...
#Sex * MAF (autosomes)
Gamma_GLM_effects_A1 <- matrix(data=NA,ncol = 3,nrow=1000)
#Sex * MAF (X)
Gamma_GLM_effects_X1 <- matrix(data=NA,ncol = 3,nrow=1000)
#Chromosome * MAF (males)
Gamma_GLM_effects_m1 <- matrix(data=NA,ncol = 3,nrow=1000)
#Chromosome * MAF (females)
Gamma_GLM_effects_f1 <- matrix(data=NA,ncol = 3,nrow=1000)
set.seed(123)
for (i in 1:1000){
  #Import Female
  fassoc1.perm <- fread(paste0(dir,"1000_permuted_phenos_weight_minus_1/",female.file.names[i]),head=T)
  fassoc1.perm$Abs_Effect <- abs(fassoc1.perm$BETA)
  fassoc1.perm$Chromosome <- fassoc1$Chromosome
  fassoc1.perm$MAF <- fassoc1$MAF
  fassoc1.perm$LD_indep <- fassoc1$LD_indep
  fassoc1.perm$Sex <- "Female"
  #Import Male
  massoc1.perm <- fread(paste0(dir,"1000_permuted_phenos_weight_minus_1/",male.file.names[i]),head=T)
  massoc1.perm$Abs_Effect <- abs(massoc1.perm$BETA)
  massoc1.perm$Chromosome <- massoc1$Chromosome
  massoc1.perm$MAF <- massoc1$MAF
  massoc1.perm$LD_indep <- massoc1$LD_indep
  massoc1.perm$Sex <- "Male"
  
  #Gamma GLM effects
  bassoc.perm <- rbind(fassoc1.perm,massoc1.perm)
  Gamma_GLM_effects_A1[i,1:3] <- as.numeric(coef(glm(data=subset(bassoc.perm,Abs_Effect>0 & LD_indep==1 & Chromosome!="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[2:4])
  Gamma_GLM_effects_X1[i,1:3] <- as.numeric(coef(glm(data=subset(bassoc.perm,Abs_Effect>0 & LD_indep==1 & Chromosome=="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[2:4])

  Gamma_GLM_effects_m1[i,1:3] <- as.numeric(coef(glm(data=subset(massoc1.perm,Abs_Effect>0 & LD_indep==1 ),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[2:4])
  Gamma_GLM_effects_f1[i,1:3] <- as.numeric(coef(glm(data=subset(fassoc1.perm,Abs_Effect>0 & LD_indep==1 ),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[2:4])
  print(i)
}

## P-values for model coefficients (autosomes)
fassoc1$Sex <- "Female"
massoc1$Sex <- "Male"
bassoc1 <- rbind(fassoc1,massoc1)
#Sex effect
sum(coef(glm(data=subset(bassoc1,Abs_Effect>0 & LD_indep==1 & Chromosome!="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[2]>Gamma_GLM_effects_A1[,1])
#p<0.001 (p=<0.001, using non-squared abs effect)
#MAF effect
sum(coef(glm(data=subset(bassoc1,Abs_Effect>0 & LD_indep==1 & Chromosome!="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[3]>Gamma_GLM_effects_A1[,2])
#p=0.760 (p=0.674, using non-squared abs effect)
#Sex*MAF effect
sum(coef(glm(data=subset(bassoc1,Abs_Effect>0 & LD_indep==1 & Chromosome!="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[4]<Gamma_GLM_effects_A1[,3])
#p=0.179 (p=0.179, using non-squared abs effect)

## P-values for model coefficients (X)
#Sex effect
sum(coef(glm(data=subset(bassoc1,Abs_Effect>0 & LD_indep==1 & Chromosome=="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[2]>Gamma_GLM_effects_X1[,1])
#p<0.001 (p=<0.001, using non-squared abs effect)
#MAF effect
sum(coef(glm(data=subset(bassoc1,Abs_Effect>0 & LD_indep==1 & Chromosome=="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[3]>Gamma_GLM_effects_X1[,2])
#p=0.191 (p=0.133, using non-squared abs effect)
#Sex*MAF effect
sum(coef(glm(data=subset(bassoc1,Abs_Effect>0 & LD_indep==1 & Chromosome=="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[4]<Gamma_GLM_effects_X1[,3])
#p=0.457 (p=0.225, using non-squared abs effect)

## ES~MAF plot (autosomes)
p1 <- ggplot(subset(bassoc1,Abs_Effect>0 & LD_indep==1& Chromosome!="X"),aes(y=Abs_Effect,x=MAF,col=Sex))+
  geom_point(alpha=0.4)+
  scale_color_manual(values=c("firebrick","dodgerblue4"))+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")
p2 <- ggplot(subset(bassoc1,Abs_Effect>0 & LD_indep==1& Chromosome!="X"),aes(y=Abs_Effect,x=MAF,col=Sex))+
  scale_color_manual(values=c("firebrick","dodgerblue4"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")+
  geom_smooth(alpha=1,method = "glm", method.args = list(family = "Gamma"))
mvsf.A <- ggdraw()+draw_plot(p1)+draw_plot(p2,x = 0.57, y = 0.57, width = 2/5, height =  2/5)

## ES~MAF plot (X chromosome)
p1 <- ggplot(subset(bassoc1,Abs_Effect>0 & LD_indep==1& Chromosome=="X"),aes(y=Abs_Effect,x=MAF,col=Sex))+
  geom_point(alpha=0.4)+
  scale_color_manual(values=c("firebrick","dodgerblue4"))+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")
p2 <- ggplot(subset(bassoc1,Abs_Effect>0 & LD_indep==1& Chromosome=="X"),aes(y=Abs_Effect,x=MAF,col=Sex))+
  scale_color_manual(values=c("firebrick","dodgerblue4"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")+
  geom_smooth(alpha=1,method = "glm", method.args = list(family = "Gamma"))
mvsf.X <- ggdraw()+draw_plot(p1)+draw_plot(p2,x = 0.57, y = 0.57, width = 2/5, height =  2/5)

## Sup Fig (scatter plot of model coefficients)

#Autosomes
Gamma_GLM_effects_A1_df <- as.data.frame(Gamma_GLM_effects_A1)
names(Gamma_GLM_effects_A1_df) <- c("Sex","MAF","Sex_by_MAF")
Gamma_GLM_effects_A1_df_melted <- melt(Gamma_GLM_effects_A1_df)
Gamma_GLM_effects_A1_df_melted$type <- "Permuted"
Gamma_GLM_effects_A1_df_melted[3001:3003,1] <- c("Sex","MAF","Sex_by_MAF")
Gamma_GLM_effects_A1_df_melted[3001:3003,2] <- coef(glm(data=subset(bassoc1,Abs_Effect>0 & LD_indep==1 & Chromosome!="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[2:4]
Gamma_GLM_effects_A1_df_melted[3001:3003,3] <- rep("True",3)
levels(Gamma_GLM_effects_A1_df_melted$variable) <- c("Sex (male)","MAF","Sex(male)-by-MAF")

#X
Gamma_GLM_effects_X1_df <- as.data.frame(Gamma_GLM_effects_X1)
names(Gamma_GLM_effects_X1_df) <- c("Sex","MAF","Sex_by_MAF")
Gamma_GLM_effects_X1_df_melted <- melt(Gamma_GLM_effects_X1_df)
Gamma_GLM_effects_X1_df_melted$type <- "Permuted"
Gamma_GLM_effects_X1_df_melted[3001:3003,1] <- c("Sex","MAF","Sex_by_MAF")
Gamma_GLM_effects_X1_df_melted[3001:3003,2] <- coef(glm(data=subset(bassoc1,Abs_Effect>0 & LD_indep==1 & Chromosome=="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[2:4]
Gamma_GLM_effects_X1_df_melted[3001:3003,3] <- rep("True",3)
levels(Gamma_GLM_effects_X1_df_melted$variable) <- c("Sex (male)","MAF","Sex(male)-by-MAF")

#Combine autosomes and X
Gamma_GLM_effects_A1_df_melted$Chromosome <- "Autosomes"
Gamma_GLM_effects_X1_df_melted$Chromosome <- "X"
Gamma_GLM_effects_df1_melted <- rbind(Gamma_GLM_effects_A1_df_melted,Gamma_GLM_effects_X1_df_melted)
pA1pX1 <- ggplot(Gamma_GLM_effects_df1_melted,aes(x=Chromosome,y=value,size=type))+
  geom_jitter(data=subset(Gamma_GLM_effects_df1_melted,type=="Permuted" & variable=="Sex(male)-by-MAF"),width=0.1,alpha=0.5,size=1.5,col="grey")+
  geom_point(data=subset(Gamma_GLM_effects_df1_melted,type=="True" & variable=="Sex(male)-by-MAF"),shape=19,size=3,col="black")+
  #geom_point(data=subset(Gamma_GLM_effects_X_df_melted,type=="True" & variable=="Sex (male)"),shape=19,size=3,col="dodgerblue4")+
  scale_size_manual(values=c(1,2))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none")+
  xlab("")+
  ylab("")

################################################
## Delta = -1####
## X vs. autosome####
################################################

#Genotypic permutations
Gamma_GLM_effects_circ_f1 <- matrix(data=NA,ncol = 3,nrow=1000)
Gamma_GLM_effects_circ_m1 <- matrix(data=NA,ncol = 3,nrow=1000)
start.point <- vector("numeric")
set.seed(123)
for (i in 1:1000){
  #Initiate X chromosome at a random starting point
  chrom.circular <- rep(fassoc1$Chromosome,2)
  start.point[i] <- sample(1:765764,1)
  fassoc1$Chromosome_perm <- chrom.circular[start.point[i]:(start.point[i]+765763)]
  massoc1$Chromosome_perm <- chrom.circular[start.point[i]:(start.point[i]+765763)]
  
  #Gamma GLM effects
  Gamma_GLM_effects_circ_f1[i,1:3] <- as.numeric(coef(glm(data=subset(fassoc1,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome_perm=="X")*MAF,family=Gamma(link="log")))[2:4])
  Gamma_GLM_effects_circ_m1[i,1:3] <- as.numeric(coef(glm(data=subset(massoc1,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome_perm=="X")*MAF,family=Gamma(link="log")))[2:4])
  print(i)
}

## P-values for each model term, autosomes vs X (females)
#Chromosome effect
Chr.effect <- as.numeric(coef(glm(data=subset(fassoc1,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[2])
sum(Chr.effect>Gamma_GLM_effects_circ_f1[,1])/1000
#p=0.33 (p=0.396 using absolute effects)
#MAF effect
MAF.effect <- as.numeric(coef(glm(data=subset(fassoc1,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[3])
sum(MAF.effect>Gamma_GLM_effects_circ_f1[,2])/1000
#p=0.859 (p=0.862 using absolute effects)
#Interaction effect
X.int.effect <- as.numeric(coef(glm(data=subset(fassoc1,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[4])
sum(X.int.effect>Gamma_GLM_effects_circ_f1[,3])/1000
#p=0.078 (p=0.099 using absolute effects)

## P-values for each model term, autosomes vs X (males)
#Chromosome effect
Chr.effect <- as.numeric(coef(glm(data=subset(massoc1,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[2])
sum(Chr.effect>Gamma_GLM_effects_circ_m1[,1])/1000
#p=0.86 (p=0.836 using absolute effects)
#MAF effect
MAF.effect <- as.numeric(coef(glm(data=subset(massoc1,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[3])
sum(MAF.effect>Gamma_GLM_effects_circ_m1[,2])/1000
#p=0.93 (p=0.925 using absolute effects)
#Interaction effect
X.int.effect <- as.numeric(coef(glm(data=subset(massoc1,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[4])
sum(X.int.effect>Gamma_GLM_effects_circ_m1[,3])/1000
#p=0.105 (p=0.147 using absolute effects)

## ES~MAF plot (females)
fassoc1$Chromosome_binary <- ifelse(fassoc1$Chromosome=='X',"X","Autosomes")
p1 <- ggplot(subset(fassoc1,Abs_Effect>0 & LD_indep==1),aes(y=Abs_Effect,x=MAF,col=Chromosome_binary))+
  geom_point(alpha=0.4)+
  scale_color_manual(values=c("forestgreen","salmon"))+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")
p2 <- ggplot(subset(fassoc1,Abs_Effect>0 & LD_indep==1),aes(y=Abs_Effect,x=MAF,col=Chromosome_binary))+
  scale_color_manual(values=c("forestgreen","salmon"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")+
  geom_smooth(alpha=1,method = "glm", method.args = list(family = "Gamma"))
XvsA.f <- ggdraw()+draw_plot(p1)+draw_plot(p2,x = 0.57, y = 0.57, width = 2/5, height =  2/5)

## ES~MAF plot (males)
massoc1$Chromosome_binary <- ifelse(massoc1$Chromosome=='X',"X","Autosomes")
p1 <- ggplot(subset(massoc1,Abs_Effect>0 & LD_indep==1),aes(y=Abs_Effect,x=MAF,col=Chromosome_binary))+
  geom_point(alpha=0.4)+
  scale_color_manual(values=c("forestgreen","salmon"))+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")
p2 <- ggplot(subset(massoc1,Abs_Effect>0 & LD_indep==1),aes(y=Abs_Effect,x=MAF,col=Chromosome_binary))+
  scale_color_manual(values=c("forestgreen","salmon"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")+
  geom_smooth(alpha=1,method = "glm", method.args = list(family = "Gamma"))
XvsA.m <- ggdraw()+draw_plot(p1)+draw_plot(p2,x = 0.57, y = 0.57, width = 2/5, height =  2/5)

## Sup Fig (scatter plot of model coefficients)

#Females
Gamma_GLM_effects_circ_f1 <- as.data.frame(Gamma_GLM_effects_circ_f1)
names(Gamma_GLM_effects_circ_f1)[1:3] <- c("Chrom.(X)","MAF","Chrom.(X)-by-MAF")
Gamma_GLM_effects_circ_f1$Mode <- "Genotype"
Gamma_GLM_effects_f1_melted1 <- melt(Gamma_GLM_effects_circ_f1)
Gamma_GLM_effects_f1 <- as.data.frame(Gamma_GLM_effects_f1)
names(Gamma_GLM_effects_f1)[1:3] <- c("Chrom.(X)","MAF","Chrom.(X)-by-MAF")
Gamma_GLM_effects_f1$Mode <- "Phenotype"
Gamma_GLM_effects_f1_melted2 <- melt(Gamma_GLM_effects_f1)
Gamma_GLM_effects_f1_melted <- rbind(Gamma_GLM_effects_f1_melted1,Gamma_GLM_effects_f1_melted2)
Gamma_GLM_effects_f1_melted$type <- "Permuted"
Gamma_GLM_effects_f1_melted[6001:6003,2] <- c("Chrom.(X)","MAF","Chrom.(X)-by-MAF")
Gamma_GLM_effects_f1_melted[6001:6003,3] <- coef(glm(data=subset(fassoc1,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[2:4]
Gamma_GLM_effects_f1_melted[6001:6003,4] <- rep("True",3)
Gamma_GLM_effects_f1_melted[6001:6003,1] <- NA
Gamma_GLM_effects_f1_melted$value <- as.numeric(as.character(Gamma_GLM_effects_f1_melted$value))

#Males
Gamma_GLM_effects_circ_m1 <- as.data.frame(Gamma_GLM_effects_circ_m1)
names(Gamma_GLM_effects_circ_m1)[1:3] <- c("Chrom.(X)","MAF","Chrom.(X)-by-MAF")
Gamma_GLM_effects_circ_m1$Mode <- "Genotype"
Gamma_GLM_effects_m1_melted1 <- melt(Gamma_GLM_effects_circ_m1)
Gamma_GLM_effects_m1 <- as.data.frame(Gamma_GLM_effects_m1)
names(Gamma_GLM_effects_m1)[1:3] <- c("Chrom.(X)","MAF","Chrom.(X)-by-MAF")
Gamma_GLM_effects_m1$Mode <- "Phenotype"
Gamma_GLM_effects_m1_melted2 <- melt(Gamma_GLM_effects_m1)
Gamma_GLM_effects_m1_melted <- rbind(Gamma_GLM_effects_m1_melted1,Gamma_GLM_effects_m1_melted2)
Gamma_GLM_effects_m1_melted$type <- "Permuted"
Gamma_GLM_effects_m1_melted[6001:6003,2] <- c("Chrom.(X)","MAF","Chrom.(X)-by-MAF")
Gamma_GLM_effects_m1_melted[6001:6003,3] <- coef(glm(data=subset(massoc1,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[2:4]
Gamma_GLM_effects_m1_melted[6001:6003,4] <- rep("True",3)
Gamma_GLM_effects_m1_melted[6001:6003,1] <- NA
Gamma_GLM_effects_m1_melted$value <- as.numeric(as.character(Gamma_GLM_effects_m1_melted$value))

#Combine females and males
Gamma_GLM_effects_f1_melted$Sex <- "Female"
Gamma_GLM_effects_m1_melted$Sex <- "Male"
Gamma_GLM_effects_fm1_melted <- rbind(Gamma_GLM_effects_f1_melted,Gamma_GLM_effects_m1_melted)
pF1pM1 <- ggplot(Gamma_GLM_effects_fm1_melted,aes(x=Sex,y=value,col=type,size=type))+
  #geom_jitter(data=subset(Gamma_GLM_effects_f_melted,Mode=="Phenotype"),width=0.1,alpha=0.5,size=1.5,col="grey")+
  geom_jitter(data=subset(Gamma_GLM_effects_fm1_melted,Mode=="Genotype" & variable=="Chrom.(X)-by-MAF"),width=0.1,alpha=0.5,size=1.5,col="grey")+
  geom_point(data=subset(Gamma_GLM_effects_fm1_melted,type=="True" & variable=="Chrom.(X)-by-MAF"),shape=19,size=3,col="black")+
  #geom_point(data=subset(Gamma_GLM_effects_fm_melted,type=="True" & variable=="Chrom.(X)"),shape=19,size=3,col="salmon")+
  scale_size_manual(values=c(1,2))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none")+
  xlab("")+
  ylab("")


################################################
## Delta = 0####
## Male vs. female####
################################################

## Import unpermuted GWAS data
fassoc0 <- read.table(paste0(dir,"female_maf0.05_weight_0.assoc"),head=T)
massoc0 <- read.table(paste0(dir,"male_maf0.05_weight_0.assoc"),head=T)
#Absolute effect size
fassoc0$Abs_Effect <- abs(fassoc0$Effect)
massoc0$Abs_Effect <- abs(massoc0$Effect)
#Chromosome coding
fassoc0$Chromosome <- as.factor(fassoc0$Chromosome)
levels(fassoc0$Chromosome) <- c("2L","2R","3L","3R","X")
massoc0$Chromosome <- as.factor(massoc0$Chromosome)
levels(massoc0$Chromosome) <- c("2L","2R","3L","3R","X")
#Title of predictor column from 'Predictor' -> 'SNP' (for LD pruning)
names(fassoc0)[2] <- "SNP"
names(massoc0)[2] <- "SNP"
#LD pruned predictors
fassoc.pruned <- read.table(paste0(dir,"female_maf0.05_pruned_r0.4_10kb.clumped"),head=T)
fassoc0$LD_indep <- ifelse(fassoc0$SNP %in% fassoc.pruned$SNP,1,0)
massoc0$LD_indep <- ifelse(massoc0$SNP %in% fassoc.pruned$SNP,1,0)

#Permuted phenotypes
female.file.names <- list.files(paste0(dir,"1000_permuted_phenos_weight_0/"),pattern = "^female_*")
male.file.names <- list.files(paste0(dir,"1000_permuted_phenos_weight_0/"),pattern = "^male_*")

#GLM coefficients for...
#Sex * MAF (autosomes)
Gamma_GLM_effects_A0 <- matrix(data=NA,ncol = 3,nrow=1000)
#Sex * MAF (X)
Gamma_GLM_effects_X0 <- matrix(data=NA,ncol = 3,nrow=1000)
#GLM coefficients for...
#Chromosome * MAF (males)
Gamma_GLM_effects_m0 <- matrix(data=NA,ncol = 3,nrow=1000)
#Chromosome * MAF (females)
Gamma_GLM_effects_f0 <- matrix(data=NA,ncol = 3,nrow=1000)
set.seed(123)
for (i in 1:1000){
  #Import Female
  fassoc0.perm <- fread(paste0(dir,"1000_permuted_phenos_weight_0/",female.file.names[i]),head=T)
  fassoc0.perm$Abs_Effect <- abs(fassoc0.perm$BETA)
  fassoc0.perm$Chromosome <- fassoc0$Chromosome
  fassoc0.perm$MAF <- fassoc0$MAF
  fassoc0.perm$LD_indep <- fassoc0$LD_indep
  fassoc0.perm$Sex <- "Female"
  #Import Male
  massoc0.perm <- fread(paste0(dir,"1000_permuted_phenos_weight_0/",male.file.names[i]),head=T)
  massoc0.perm$Abs_Effect <- abs(massoc0.perm$BETA)
  massoc0.perm$Chromosome <- massoc0$Chromosome
  massoc0.perm$MAF <- massoc0$MAF
  massoc0.perm$LD_indep <- massoc0$LD_indep
  massoc0.perm$Sex <- "Male"
  
  #Gamma GLM effects
  bassoc.perm <- rbind(fassoc0.perm,massoc0.perm)
  Gamma_GLM_effects_A0[i,1:3] <- as.numeric(coef(glm(data=subset(bassoc.perm,Abs_Effect>0 & LD_indep==1 & Chromosome!="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[2:4])
  Gamma_GLM_effects_X0[i,1:3] <- as.numeric(coef(glm(data=subset(bassoc.perm,Abs_Effect>0 & LD_indep==1 & Chromosome=="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[2:4])

  Gamma_GLM_effects_m0[i,1:3] <- as.numeric(coef(glm(data=subset(massoc0.perm,Abs_Effect>0 & LD_indep==1 ),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[2:4])
  Gamma_GLM_effects_f0[i,1:3] <- as.numeric(coef(glm(data=subset(fassoc0.perm,Abs_Effect>0 & LD_indep==1 ),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[2:4])
  print(i)
}

## P-values for model coefficients, male vs. female (autosomes)
fassoc0$Sex <- "Female"
massoc0$Sex <- "Male"
bassoc0 <- rbind(fassoc0,massoc0)
#Sex effect
sum(coef(glm(data=subset(bassoc0,Abs_Effect>0 & LD_indep==1 & Chromosome!="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[2]>Gamma_GLM_effects_A0[,1])
#p<0.001 (p<0.001 using absolute effects)
#MAF effect
sum(coef(glm(data=subset(bassoc0,Abs_Effect>0 & LD_indep==1 & Chromosome!="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[3]>Gamma_GLM_effects_A0[,2])
#p=0.829 (p=0.737 using absolute effects)
#Sex*MAF effect
sum(coef(glm(data=subset(bassoc0,Abs_Effect>0 & LD_indep==1 & Chromosome!="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[4]<Gamma_GLM_effects_A0[,3])
#p=0.285 (p=0.201 using absolute effects)

## P-values for model coefficients, male vs. female (X)
#Sex effect
sum(coef(glm(data=subset(bassoc0,Abs_Effect>0 & LD_indep==1 & Chromosome=="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[2]>Gamma_GLM_effects_X0[,1])
#p<0.001 (p<0.001 using absolute effects)
#MAF effect
sum(coef(glm(data=subset(bassoc0,Abs_Effect>0 & LD_indep==1 & Chromosome=="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[3]>Gamma_GLM_effects_X0[,2])
#p=0.183 (p=0.120 using absolute effects)
#Sex*MAF effect
sum(coef(glm(data=subset(bassoc0,Abs_Effect>0 & LD_indep==1 & Chromosome=="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[4]<Gamma_GLM_effects_X0[,3])
#p=0.430 (p=0.211 using absolute effects)

## ES~MAF plot (autosomes)
p1 <- ggplot(subset(bassoc0,Abs_Effect>0 & LD_indep==1& Chromosome!="X"),aes(y=Abs_Effect,x=MAF,col=Sex))+
  geom_point(alpha=0.4)+
  scale_color_manual(values=c("firebrick","dodgerblue4"))+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")
p2 <- ggplot(subset(bassoc0,Abs_Effect>0 & LD_indep==1& Chromosome!="X"),aes(y=Abs_Effect,x=MAF,col=Sex))+
  scale_color_manual(values=c("firebrick","dodgerblue4"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")+
  geom_smooth(alpha=1,method = "glm", method.args = list(family = "Gamma"))
mvsf.A <- ggdraw()+draw_plot(p1)+draw_plot(p2,x = 0.57, y = 0.57, width = 2/5, height =  2/5)

## ES~MAF plot (X chromosome)
p1 <- ggplot(subset(bassoc0,Abs_Effect>0 & LD_indep==1& Chromosome=="X"),aes(y=Abs_Effect,x=MAF,col=Sex))+
  geom_point(alpha=0.4)+
  scale_color_manual(values=c("firebrick","dodgerblue4"))+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")
p2 <- ggplot(subset(bassoc0,Abs_Effect>0 & LD_indep==1& Chromosome=="X"),aes(y=Abs_Effect,x=MAF,col=Sex))+
  scale_color_manual(values=c("firebrick","dodgerblue4"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")+
  geom_smooth(alpha=1,method = "glm", method.args = list(family = "Gamma"))
mvsf.X <- ggdraw()+draw_plot(p1)+draw_plot(p2,x = 0.57, y = 0.57, width = 2/5, height =  2/5)

## Sup Fig (scatter plot of model coefficients)

#Autosomes
Gamma_GLM_effects_A0_df <- as.data.frame(Gamma_GLM_effects_A0)
names(Gamma_GLM_effects_A0_df) <- c("Sex","MAF","Sex_by_MAF")
Gamma_GLM_effects_A0_df_melted <- melt(Gamma_GLM_effects_A0_df)
Gamma_GLM_effects_A0_df_melted$type <- "Permuted"
Gamma_GLM_effects_A0_df_melted[3001:3003,1] <- c("Sex","MAF","Sex_by_MAF")
Gamma_GLM_effects_A0_df_melted[3001:3003,2] <- coef(glm(data=subset(bassoc0,Abs_Effect>0 & LD_indep==1 & Chromosome!="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[2:4]
Gamma_GLM_effects_A0_df_melted[3001:3003,3] <- rep("True",3)
levels(Gamma_GLM_effects_A0_df_melted$variable) <- c("Sex (male)","MAF","Sex(male)-by-MAF")

#X
Gamma_GLM_effects_X0_df <- as.data.frame(Gamma_GLM_effects_X0)
names(Gamma_GLM_effects_X0_df) <- c("Sex","MAF","Sex_by_MAF")
Gamma_GLM_effects_X0_df_melted <- melt(Gamma_GLM_effects_X0_df)
Gamma_GLM_effects_X0_df_melted$type <- "Permuted"
Gamma_GLM_effects_X0_df_melted[3001:3003,1] <- c("Sex","MAF","Sex_by_MAF")
Gamma_GLM_effects_X0_df_melted[3001:3003,2] <- coef(glm(data=subset(bassoc0,Abs_Effect>0 & LD_indep==1 & Chromosome=="X"),Abs_Effect~Sex*MAF,family=Gamma(link="log")))[2:4]
Gamma_GLM_effects_X0_df_melted[3001:3003,3] <- rep("True",3)
levels(Gamma_GLM_effects_X0_df_melted$variable) <- c("Sex (male)","MAF","Sex(male)-by-MAF")

#Combine autosomes and X
Gamma_GLM_effects_A0_df_melted$Chromosome <- "Autosomes"
Gamma_GLM_effects_X0_df_melted$Chromosome <- "X"
Gamma_GLM_effects_df0_melted <- rbind(Gamma_GLM_effects_A0_df_melted,Gamma_GLM_effects_X0_df_melted)
pA0pX0 <- ggplot(Gamma_GLM_effects_df0_melted,aes(x=Chromosome,y=value,size=type))+
  geom_jitter(data=subset(Gamma_GLM_effects_df0_melted,type=="Permuted" & variable=="Sex(male)-by-MAF"),width=0.1,alpha=0.5,size=1.5,col="grey")+
  geom_point(data=subset(Gamma_GLM_effects_df0_melted,type=="True"& variable=="Sex(male)-by-MAF"),shape=19,size=3,col="black")+
  #geom_point(data=subset(Gamma_GLM_effects_X_df_melted,type=="True" & variable=="Sex (male)"),shape=19,size=3,col="dodgerblue4")+
  scale_size_manual(values=c(1,2))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none")+
  xlab("")+
  ylab("")


################################################
## Delta = 0####
## X vs. autosome####
################################################

#Genotypic permutations
Gamma_GLM_effects_circ_f0 <- matrix(data=NA,ncol = 3,nrow=1000)
Gamma_GLM_effects_circ_m0 <- matrix(data=NA,ncol = 3,nrow=1000)
start.point <- vector("numeric")
set.seed(123)
for (i in 1:1000){
  #Initiate X chromosome at a random starting point
  chrom.circular <- rep(fassoc0$Chromosome,2)
  start.point[i] <- sample(1:765764,1)
  fassoc0$Chromosome_perm <- chrom.circular[start.point[i]:(start.point[i]+765763)]
  massoc0$Chromosome_perm <- chrom.circular[start.point[i]:(start.point[i]+765763)]
  
  #Gamma GLM effects
  Gamma_GLM_effects_circ_f0[i,1:3] <- as.numeric(coef(glm(data=subset(fassoc0,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome_perm=="X")*MAF,family=Gamma(link="log")))[2:4])
  Gamma_GLM_effects_circ_m0[i,1:3] <- as.numeric(coef(glm(data=subset(massoc0,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome_perm=="X")*MAF,family=Gamma(link="log")))[2:4])
  
  print(i)
}

## P-values for each model term, autosomes vs X (females)
#Chromosome effect
Chr.effect <- as.numeric(coef(glm(data=subset(fassoc0,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[2])
sum(Chr.effect>Gamma_GLM_effects_circ_f0[,1])/1000
#p=0.406 (p=0.478 using absolute effects)
#MAF effect
MAF.effect <- as.numeric(coef(glm(data=subset(fassoc0,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[3])
sum(MAF.effect>Gamma_GLM_effects_circ_f0[,2])/1000
#p=0.936 (p=0.934 using absolute effects)
#Interaction effect
X.int.effect <- as.numeric(coef(glm(data=subset(fassoc0,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[4])
sum(X.int.effect>Gamma_GLM_effects_circ_f0[,3])/1000
#p=0.032 (p=0.048 using absolute effects)

## P-values for each model term, autosomes vs X (males)
#Chromosome effect
Chr.effect <- as.numeric(coef(glm(data=subset(massoc0,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[2])
sum(Chr.effect>Gamma_GLM_effects_circ_m0[,1])/1000
#p=0.854 (p=0.831 using absolute effects)
#MAF effect
MAF.effect <- as.numeric(coef(glm(data=subset(massoc0,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[3])
sum(MAF.effect>Gamma_GLM_effects_circ_m0[,2])/1000
#p=0.933 (p=0.925 using absolute effects)
#Interaction effect
X.int.effect <- as.numeric(coef(glm(data=subset(massoc0,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[4])
sum(X.int.effect>Gamma_GLM_effects_circ_m0[,3])/1000
#p=0.094 (p=0.139 using absolute effects)

## ES~MAF plot (females)
fassoc0$Chromosome_binary <- ifelse(fassoc0$Chromosome=='X',"X","Autosomes")
p1 <- ggplot(subset(fassoc0,Abs_Effect>0 & LD_indep==1),aes(y=Abs_Effect,x=MAF,col=Chromosome_binary))+
  geom_point(alpha=0.4)+
  scale_color_manual(values=c("forestgreen","salmon"))+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")
p2 <- ggplot(subset(fassoc0,Abs_Effect>0 & LD_indep==1),aes(y=Abs_Effect,x=MAF,col=Chromosome_binary))+
  scale_color_manual(values=c("forestgreen","salmon"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")+
  geom_smooth(alpha=1,method = "glm", method.args = list(family = "Gamma"))
XvsA.f <- ggdraw()+draw_plot(p1)+draw_plot(p2,x = 0.57, y = 0.57, width = 2/5, height =  2/5)

## ES~MAF plot (males)
massoc0$Chromosome_binary <- ifelse(massoc0$Chromosome=='X',"X","Autosomes")
p1 <- ggplot(subset(massoc0,Abs_Effect>0 & LD_indep==1),aes(y=Abs_Effect,x=MAF,col=Chromosome_binary))+
  geom_point(alpha=0.4)+
  scale_color_manual(values=c("forestgreen","salmon"))+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")
p2 <- ggplot(subset(massoc0,Abs_Effect>0 & LD_indep==1),aes(y=Abs_Effect,x=MAF,col=Chromosome_binary))+
  scale_color_manual(values=c("forestgreen","salmon"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none")+
  ylab(expression(paste("|",alpha,"|")))+
  xlab("MAF")+
  geom_smooth(alpha=1,method = "glm", method.args = list(family = "Gamma"))
XvsA.m <- ggdraw()+draw_plot(p1)+draw_plot(p2,x = 0.57, y = 0.57, width = 2/5, height =  2/5)

## Sup Fig (scatter plot of model coefficients)

#Females
Gamma_GLM_effects_circ_f0 <- as.data.frame(Gamma_GLM_effects_circ_f0)
names(Gamma_GLM_effects_circ_f0)[1:3] <- c("Chrom.(X)","MAF","Chrom.(X)-by-MAF")
Gamma_GLM_effects_circ_f0$Mode <- "Genotype"
Gamma_GLM_effects_f0_melted1 <- melt(Gamma_GLM_effects_circ_f0)
Gamma_GLM_effects_f0 <- as.data.frame(Gamma_GLM_effects_f0)
names(Gamma_GLM_effects_f0)[1:3] <- c("Chrom.(X)","MAF","Chrom.(X)-by-MAF")
Gamma_GLM_effects_f0$Mode <- "Phenotype"
Gamma_GLM_effects_f0_melted2 <- melt(Gamma_GLM_effects_f0)
Gamma_GLM_effects_f0_melted <- rbind(Gamma_GLM_effects_f0_melted1,Gamma_GLM_effects_f0_melted2)
Gamma_GLM_effects_f0_melted$type <- "Permuted"
Gamma_GLM_effects_f0_melted[6001:6003,2] <- c("Chrom.(X)","MAF","Chrom.(X)-by-MAF")
Gamma_GLM_effects_f0_melted[6001:6003,3] <- coef(glm(data=subset(fassoc0,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[2:4]
Gamma_GLM_effects_f0_melted[6001:6003,4] <- rep("True",3)
Gamma_GLM_effects_f0_melted[6001:6003,1] <- NA
Gamma_GLM_effects_f0_melted$value <- as.numeric(as.character(Gamma_GLM_effects_f0_melted$value))

#Males
Gamma_GLM_effects_circ_m0 <- as.data.frame(Gamma_GLM_effects_circ_m0)
names(Gamma_GLM_effects_circ_m0) <- c("Chromosome_region","MAF","Region_by_MAF")
Gamma_GLM_effects_circ_m0 <- as.data.frame(Gamma_GLM_effects_circ_m0)
names(Gamma_GLM_effects_circ_m0)[1:3] <- c("Chrom.(X)","MAF","Chrom.(X)-by-MAF")
Gamma_GLM_effects_circ_m0$Mode <- "Genotype"
Gamma_GLM_effects_m0_melted1 <- melt(Gamma_GLM_effects_circ_m0)
Gamma_GLM_effects_m0 <- as.data.frame(Gamma_GLM_effects_m0)
names(Gamma_GLM_effects_m0)[1:3] <- c("Chrom.(X)","MAF","Chrom.(X)-by-MAF")
Gamma_GLM_effects_m0$Mode <- "Phenotype"
Gamma_GLM_effects_m0_melted2 <- melt(Gamma_GLM_effects_m0)
Gamma_GLM_effects_m0_melted <- rbind(Gamma_GLM_effects_m0_melted1,Gamma_GLM_effects_m0_melted2)
Gamma_GLM_effects_m0_melted$type <- "Permuted"
Gamma_GLM_effects_m0_melted[6001:6003,2] <- c("Chrom.(X)","MAF","Chrom.(X)-by-MAF")
Gamma_GLM_effects_m0_melted[6001:6003,3] <- coef(glm(data=subset(massoc0,Abs_Effect>0 & LD_indep==1),Abs_Effect~factor(Chromosome=="X")*MAF,family=Gamma(link="log")))[2:4]
Gamma_GLM_effects_m0_melted[6001:6003,4] <- rep("True",3)
Gamma_GLM_effects_m0_melted[6001:6003,1] <- NA
Gamma_GLM_effects_m0_melted$value <- as.numeric(as.character(Gamma_GLM_effects_m0_melted$value))

#Combine males and females
Gamma_GLM_effects_f0_melted$Sex <- "Female"
Gamma_GLM_effects_m0_melted$Sex <- "Male"
Gamma_GLM_effects_fm0_melted <- rbind(Gamma_GLM_effects_f0_melted,Gamma_GLM_effects_m0_melted)
pF0pM0 <- ggplot(Gamma_GLM_effects_fm0_melted,aes(x=Sex,y=value,col=type,size=type))+
  #geom_jitter(data=subset(Gamma_GLM_effects_f_melted,Mode=="Phenotype"),width=0.1,alpha=0.5,size=1.5,col="grey")+
  geom_jitter(data=subset(Gamma_GLM_effects_fm0_melted,Mode=="Genotype" & variable=="Chrom.(X)-by-MAF"),width=0.1,alpha=0.5,size=1.5,col="grey")+
  geom_point(data=subset(Gamma_GLM_effects_fm0_melted,type=="True" & variable=="Chrom.(X)-by-MAF"),shape=19,size=3,col="black")+
  #geom_point(data=subset(Gamma_GLM_effects_fm_melted,type=="True" & variable=="Chrom.(X)"),shape=19,size=3,col="salmon")+
  scale_size_manual(values=c(1,2))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none")+
  xlab("")+
  ylab("")

################################################
## Combined deltas
################################################

## Combine SupFigs across all alphas
plot_grid(pApX,pA1pX1,pA0pX0,labels=c("Delta = -0.25","Delta = -1","Delta = 0"),scale = 0.9,ncol=3)
plot_grid(pFpM,pF1pM1,pF0pM0,labels=c("Delta = -0.25","Delta = -1","Delta = 0"),scale = 0.9,ncol=3)
