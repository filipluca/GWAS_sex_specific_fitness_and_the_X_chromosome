library(ggplot2)
library(cowplot)

dir <- "~/Dropbox/dropbox_work/data/GWAS_m_f_fitness/"

############################################################
## Correlation between gene length and gene-based heritability####
############################################################

## Females
gbat.f <- read.table(paste0(dir,"gbat_female_maf0.05/remls.1"),head=T)
ensgenes <- read.table(paste0(dir,"ensgenes.txt"))
gbat.f.X <- subset(gbat.f,Gene_Name %in% subset(ensgenes,V2==5)$V1)
gbat.f.A <- subset(gbat.f,Gene_Name %in% subset(ensgenes,V2!=5)$V1)
gbat.f$Chromosome <- ifelse(gbat.f$Gene_Name %in% subset(ensgenes,V2==5)$V1,"X","Autosomes")

## Whole-genome
#Permuted correlation
cors.f <- vector("numeric")
for (i in 12:1011){
  cors.f[i] <- as.numeric(cor.test(gbat.f$Length,gbat.f[,i],method="spearman")$estimate)
}
cors.f <- as.data.frame(cors.f)
names(cors.f) <- "Correlations"
cors.f$Type <- "Permuted"

#Observed correlation
t.cors.f <- as.numeric(cor.test(gbat.f$Length,gbat.f$Heritability,method="spearman")$estimate)
cors.f[1012,1] <- t.cors.f
cors.f[1012,2] <- "True"
cors.f$Sex <- "Female"
#p-value
sum(t.cors.f<cors.f$Correlations,na.rm=T)/1000

## Autosomes
#Permuted correlation
cors.f.A <- vector("numeric")
for (i in 12:1011){
  cors.f.A[i] <- as.numeric(cor.test(gbat.f.A$Length,gbat.f.A[,i],method="spearman")$estimate)
}
cors.f.A <- as.data.frame(cors.f.A)
names(cors.f.A) <- "Correlations"
cors.f.A$Type <- "Permuted"
#Observed correlation
t.cors.f.A <- as.numeric(cor.test(gbat.f.A$Length,gbat.f.A$Heritability,method="spearman")$estimate)
cors.f.A[1012,1] <- t.cors.f.A
cors.f.A[1012,2] <- "True"
cors.f.A$Sex <- "Female"
cors.f.A$Chromosome <- "Autosomes"
#p-value
sum(t.cors.f.A<cors.f.A$Correlations,na.rm=T)/1000
#0.000

#X chromosome
cors.f.X <- vector("numeric")
for (i in 12:1011){
  cors.f.X[i] <- as.numeric(cor.test(gbat.f.X$Length,gbat.f.X[,i],method="spearman")$estimate)
}
cors.f.X <- as.data.frame(cors.f.X)
names(cors.f.X) <- "Correlations"
cors.f.X$Type <- "Permuted"
#Observed correlation
t.cors.f.X <- as.numeric(cor.test(gbat.f.X$Length,gbat.f.X$Heritability,method="spearman")$estimate)
cors.f.X[1012,1] <- t.cors.f.X
cors.f.X[1012,2] <- "True"
cors.f.X$Sex <- "Female"
cors.f.X$Chromosome <- "X"
#p-value
sum(t.cors.f.X<cors.f.X$Correlations,na.rm=T)/1000
#0.024


##Males
gbat.m <- read.table(paste0(dir,"gbat_male_maf0.05/remls.1"),head=T)
ensgenes <- read.table(paste0(dir,"ensgenes.txt"))
gbat.m.X <- subset(gbat.m,Gene_Name %in% subset(ensgenes,V2==5)$V1)
gbat.m.A <- subset(gbat.m,Gene_Name %in% subset(ensgenes,V2!=5)$V1)
gbat.m$Chromosome <- ifelse(gbat.m$Gene_Name %in% subset(ensgenes,V2==5)$V1,"X","Autosomes")

## Whole-genome
#Permuted correlation
cors.m <- vector("numeric")
for (i in 12:1011){
  cors.m[i] <- as.numeric(cor.test(gbat.m$Length,gbat.m[,i],method="spearman")$estimate)
}
cors.m <- as.data.frame(cors.m)
names(cors.m) <- "Correlations"
cors.m$Type <- "Permuted"
#Observed correlation
t.cors.m <- as.numeric(cor.test(gbat.m$Length,gbat.m$Heritability,method="spearman")$estimate)
cors.m[1012,1] <- t.cors.m
cors.m[1012,2] <- "True"
cors.m$Sex <- "Male"
#p-value
sum(t.cors.m<cors.m$Correlations,na.rm=T)/1000
#0.000

## Autosomes
#Permuted correlation
cors.m.A <- vector("numeric")
for (i in 12:1011){
  cors.m.A[i] <- as.numeric(cor.test(gbat.m.A$Length,gbat.m.A[,i],method="spearman")$estimate)
}
cors.m.A <- as.data.frame(cors.m.A)
names(cors.m.A) <- "Correlations"
cors.m.A$Type <- "Permuted"
#Observed correlation
t.cors.m.A <- as.numeric(cor.test(gbat.m.A$Length,gbat.m.A$Heritability,method="spearman")$estimate)
cors.m.A[1012,1] <- t.cors.m.A
cors.m.A[1012,2] <- "True"
cors.m.A$Sex <- "Male"
cors.m.A$Chromosome <- "Autosomes"
#p-value
sum(t.cors.m.A<cors.m.A$Correlations,na.rm=T)/1000
#0.000

## X
#Permuted correlation
cors.m.X <- vector("numeric")
for (i in 12:1011){
  cors.m.X[i] <- as.numeric(cor.test(gbat.m.X$Length,gbat.m.X[,i],method="spearman")$estimate)
}
cors.m.X <- as.data.frame(cors.m.X)
names(cors.m.X) <- "Correlations"
cors.m.X$Type <- "Permuted"
#Observed correlation
t.cors.m.X <- as.numeric(cor.test(gbat.m.X$Length,gbat.m.X$Heritability,method="spearman")$estimate)
cors.m.X[1012,1] <- t.cors.m.X
cors.m.X[1012,2] <- "True"
cors.m.X$Sex <- "Male"
cors.m.X$Chromosome <- "X"
#p-value
sum(t.cors.m.X<cors.m.X$Correlations,na.rm=T)/1000
#0.849


#################################################################################
## Prepare chunks ####
#################################################################################

fassoc <- read.table(paste0(dir,"female_maf0.05.assoc"),h=T)
fassoc <- fassoc[,c("Predictor","Chromosome","Basepair")]

rownames_2L <- as.numeric(row.names(subset(fassoc,Chromosome==1)))
rownames_2R <- as.numeric(row.names(subset(fassoc,Chromosome==2)))
rownames_3L <- as.numeric(row.names(subset(fassoc,Chromosome==3)))
rownames_3R <- as.numeric(row.names(subset(fassoc,Chromosome==4)))
rownames_X <- as.numeric(row.names(subset(fassoc,Chromosome==5)))

set.seed(123)
for (i in 1:1000){
  random_breakpoint_2L <- sort(sample(rownames_2L,replace=F,size = 499))
  random_breakpoint_2R <- sort(sample(rownames_2R,replace=F,size = 499))
  random_breakpoint_3L <- sort(sample(rownames_3L,replace=F,size = 499))
  random_breakpoint_3R <- sort(sample(rownames_3R,replace=F,size = 499))
  random_breakpoint_X <- sort(sample(rownames_X,replace=F,size = 499))
  
  #Make chunks file for autosomes
  chunks1_2L <- rbind(fassoc[head(rownames_2L,1),],fassoc[random_breakpoint_2L,],fassoc[tail(rownames_2L,1),])
  chunks2_2L <- rbind(fassoc[random_breakpoint_2L,],fassoc[tail(rownames_2L,1),])
  chunks_2L <- cbind(chunks1_2L[1:500,],chunks2_2L[c("Basepair")])
  
  chunks1_2R <- rbind(fassoc[head(rownames_2R,1),],fassoc[random_breakpoint_2R,],fassoc[tail(rownames_2R,1),])
  chunks2_2R <- rbind(fassoc[random_breakpoint_2R,],fassoc[tail(rownames_2R,1),])
  chunks_2R <- cbind(chunks1_2R[1:500,],chunks2_2R[c("Basepair")])
  
  chunks1_3L <- rbind(fassoc[head(rownames_3L,1),],fassoc[random_breakpoint_3L,],fassoc[tail(rownames_3L,1),])
  chunks2_3L <- rbind(fassoc[random_breakpoint_3L,],fassoc[tail(rownames_3L,1),])
  chunks_3L <- cbind(chunks1_3L[1:500,],chunks2_3L[c("Basepair")])
  
  chunks1_3R <- rbind(fassoc[head(rownames_3R,1),],fassoc[random_breakpoint_3R,],fassoc[tail(rownames_3R,1),])
  chunks2_3R <- rbind(fassoc[random_breakpoint_3R,],fassoc[tail(rownames_3R,1),])
  chunks_3R <- cbind(chunks1_3R[1:500,],chunks2_3R[c("Basepair")])
  
  chunks1_X <- rbind(fassoc[head(rownames_X,1),],fassoc[random_breakpoint_X,],fassoc[tail(rownames_X,1),])
  chunks2_X <- rbind(fassoc[random_breakpoint_X,],fassoc[tail(rownames_X,1),])
  chunks_X <- cbind(chunks1_X[1:500,],chunks2_X[c("Basepair")])
  
  chunks_all <- rbind(chunks_2L,chunks_2R,chunks_3L,chunks_3R,chunks_X)
  print(i)
  write.table(chunks_all,paste0(dir,"Dmel_chr_chunks/Dmel_chr_chunks",i,".txt"),quote=F,col.names=F,row.names=F)
}


############################################################
## Correlation between chunk length and chunk-based heritability####
############################################################

#Males
cors.m.A <- vector("numeric")
cors.m.X <- vector("numeric")
cv.m.A <- vector("numeric")
cv.m.X <- vector("numeric")

for (i in 1:1000){
  cbat.m <- read.table(paste0(dir,"chunk_bat_male_maf0.05/remls",i),h=T)
  cbat.m.A <- cbat.m[1:2000,]
  cbat.m.X <- cbat.m[2001:2500,]
  cors.m.A[i] <- cor.test(cbat.m.A$Length,cbat.m.A$Heritability,method="spearman")$estimate
  cors.m.X[i] <- cor.test(cbat.m.X$Length,cbat.m.X$Heritability,method="spearman")$estimate
  cv.m.A[i] <- sqrt(var(cbat.m.A$Heritability))/mean(cbat.m.A$Heritability)
  cv.m.X[i] <- sqrt(var(cbat.m.X$Heritability))/mean(cbat.m.X$Heritability)
}

#Females
cors.f.A <- vector("numeric")
cors.f.X <- vector("numeric")
cv.f.A <- vector("numeric")
cv.f.X <- vector("numeric")

for (i in 1:1000){
  cbat.f <- read.table(paste0(dir,"chunk_bat_female_maf0.05/remls",i),h=T)
  cbat.f.A <- cbat.f[1:2000,]
  cbat.f.X <- cbat.f[2001:2500,]
  cors.f.A[i] <- cor.test(cbat.f.A$Length,cbat.f.A$Heritability,method="spearman")$estimate
  cors.f.X[i] <- cor.test(cbat.f.X$Length,cbat.f.X$Heritability,method="spearman")$estimate
  cv.f.A[i] <- sqrt(var(cbat.f.A$Heritability))/mean(cbat.f.A$Heritability)
  cv.f.X[i] <- sqrt(var(cbat.f.X$Heritability))/mean(cbat.f.X$Heritability)
}

## Plotting chunks
stats.m.X <- as.data.frame(cbind(cors.m.X,cv.m.X))
stats.m.A <- as.data.frame(cbind(cors.m.A,cv.m.A))
names(stats.m.X) <- c("Correlation","CV")
names(stats.m.A) <- c("Correlation","CV")
stats.m.X$Sex <- "Male"
stats.m.A$Sex <- "Male"
stats.m.X$Chromosome <- "X"
stats.m.A$Chromosome <- "Autosomes"

stats.f.X <- as.data.frame(cbind(cors.f.X,cv.f.X))
stats.f.A <- as.data.frame(cbind(cors.f.A,cv.f.A))
names(stats.f.X) <- c("Correlation","CV")
names(stats.f.A) <- c("Correlation","CV")
stats.f.X$Sex <- "Female"
stats.f.A$Sex <- "Female"
stats.f.X$Chromosome <- "X"
stats.f.A$Chromosome <- "Autosomes"

#Difference between male and female correlations
stats.diff.A <- stats.f.A$Correlation-stats.m.A$Correlation
c(median(stats.diff.A),quantile(stats.diff.A,0.025),quantile(stats.diff.A,0.975))

stats.diff.X <- stats.f.X$Correlation-stats.m.X$Correlation
c(median(stats.diff.X),quantile(stats.diff.X,0.025),quantile(stats.diff.X,0.975))

#p-values
#Females,A
sum(stats.f.A$Correlation<0)/1000
median(stats.f.A$Correlation)
quantile(stats.f.A$Correlation,0.975)
quantile(stats.f.A$Correlation,0.025)
#Females,X
sum(stats.f.X$Correlation<0)/1000
median(stats.f.X$Correlation)
quantile(stats.f.X$Correlation,0.975)
quantile(stats.f.X$Correlation,0.025)

#Males,A
sum(stats.m.A$Correlation<0)/1000
median(stats.m.A$Correlation)
quantile(stats.m.A$Correlation,0.975)
quantile(stats.m.A$Correlation,0.025)
#Males,X
sum(stats.m.X$Correlation<0)/1000
median(stats.m.X$Correlation)
quantile(stats.m.X$Correlation,0.975)
quantile(stats.m.X$Correlation,0.025)


##########
#Plots ####
##########

f1 <- ggplot(gbat.f,aes(log10(Length),y=log10(Heritability)))+
  geom_point(size=1,alpha=0.3,col="firebrick4")+
  theme_bw()+
  stat_smooth(method="lm",col="black")+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),strip.background = element_rect(fill="white"),strip.text = element_text(size=15))+
  ylab(expression(paste(log[10],"(",h[SNP]^2,")")))+
  xlab(expression(paste(log[10],"(Gene length)")))+
  facet_wrap(~Chromosome)

m1 <- ggplot(gbat.m,aes(log10(Length),y=log10(Heritability)))+
  geom_point(size=1,alpha=0.3,col="dodgerblue4")+
  theme_bw()+
  stat_smooth(method="lm",col="black")+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),strip.background = element_rect(fill="white"),strip.text = element_text(size=15))+
  ylab(expression(paste(log[10],"(",h[SNP]^2,")")))+
  xlab(expression(paste(log[10],"(Gene length)")))+
  facet_wrap(~Chromosome)

cors.f.b <- rbind(cors.f.A,cors.f.X)
f2 <- ggplot(cors.f.b,aes(x=Chromosome,y=Correlations,col=Type,size=Type))+
  geom_jitter(data=subset(cors.f.b,Type=="Permuted" & Sex=="Female"),width=0.1,alpha=0.5,size=1.5,col="grey")+
  #geom_point(data=subset(cors.b,Type=="True" & Sex=="Male"),shape=19,size=3,col="dodgerblue4")+
  geom_point(data=subset(cors.f.b,Type=="True"& Sex=="Female"),shape=19,size=3,col="firebrick3")+
  scale_size_manual(values=c(1,2))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none",strip.background = element_rect(fill="white"),strip.text = element_text(size=15))+
  ylab(expression(paste(rho,"(Gene length, h"[SNP]^2,")")))+
  xlab("")+
  ylim(c(-0.05,0.11))


cors.m.b <- rbind(cors.m.A,cors.m.X)
m2 <- ggplot(cors.m.b,aes(x=Chromosome,y=Correlations,col=Type,size=Type))+
  geom_jitter(data=subset(cors.m.b,Type=="Permuted"  & Sex=="Male"),width=0.1,alpha=0.5,size=1.5,col="grey")+
  geom_point(data=subset(cors.m.b,Type=="True" & Sex=="Male"),shape=19,size=3,col="dodgerblue4")+
  #geom_point(data=subset(cors.b,Type=="True"& Sex=="Female"),shape=19,size=3,col="firebrick3")+
  scale_size_manual(values=c(1,2))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none",strip.background = element_rect(fill="white"),strip.text = element_text(size=15))+
  ylab(expression(paste(rho,"(Gene length, h"[SNP]^2,")")))+
  xlab("")+
  ylim(c(-0.05,0.11))

stats.b <- rbind(stats.m.A,stats.m.X,stats.f.A,stats.f.X)

ggplot(stats.b,aes(x=Chromosome,y=Correlation,col=Sex))+
  stat_summary(fun.ymin = function(x) quantile(x,0.025), fun.ymax = function(x) quantile(x,0.975),width=0.3,geom="errorbar",size=0.8)+
  stat_summary(fun.y=mean, geom="point",size=3)+
  #geom_jitter(alpha=0.15,width=0.2)+
  scale_color_manual(values=c("dodgerblue4","firebrick"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none",strip.background = element_rect(fill="white"),strip.text = element_text(size=15))+
  ylab(expression(paste(rho,"(Partition length, h"[SNP]^2,")")))+
  xlab("")+
  facet_wrap(~Sex)+
  geom_hline(yintercept = 0,linetype="dashed")

m3 <- ggplot(subset(stats.b,Sex=="Male"),aes(x=Chromosome,y=Correlation,col=Sex))+
  geom_jitter(alpha=0.3,width=0.2,col="dodgerblue4")+
  stat_summary(fun.ymin = function(x) quantile(x,0.025), fun.ymax = function(x) quantile(x,0.975),width=0.3,geom="errorbar",size=0.8)+
  stat_summary(fun.y=mean, geom="point",size=3)+
  scale_color_manual(values=c("black"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none",strip.background = element_rect(fill="white"),strip.text = element_text(size=15))+
  ylab(expression(paste(rho,"(Partition length, h"[SNP]^2,")")))+
  xlab("")+
  #facet_wrap(~Sex)+
  geom_hline(yintercept = 0,linetype="dashed")

f3 <- ggplot(subset(stats.b,Sex=="Female"),aes(x=Chromosome,y=Correlation,col=Sex))+
  geom_jitter(alpha=0.3,width=0.2,col="firebrick")+
  stat_summary(fun.ymin = function(x) quantile(x,0.025), fun.ymax = function(x) quantile(x,0.975),width=0.3,geom="errorbar",size=0.8)+
  stat_summary(fun.y=mean, geom="point",size=3)+
  scale_color_manual(values=c("black"))+
  theme_bw()+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none",strip.background = element_rect(fill="white"),strip.text = element_text(size=15))+
  ylab(expression(paste(rho,"(Partition length, h"[SNP]^2,")")))+
  xlab("")+
  #facet_wrap(~Sex)+
  geom_hline(yintercept = 0,linetype="dashed")


#All plots together
#plot_grid(f1,f2,m1,m2,ncol=2,scale = 0.9,labels = c("Female",NA,"Male",NA),rel_widths = c(2,1,2,1),hjust = 0)
#plot_grid(f1,m1,ncol=1,scale = 0.9,labels = c("Female","Male"),rel_widths = c(2,1,2,1),hjust = 0)
plot_grid(f1,f3,m1,m3,ncol=2,scale = 0.9,labels = c("Female",NA,"Male",NA),rel_widths = c(2,1,2,1),hjust = 0)



