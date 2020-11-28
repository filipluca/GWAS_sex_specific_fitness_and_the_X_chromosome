library(CMplot)
library(qqman)
library(Hmisc)
library(reshape2)
library(ggplot2)
library(data.table)
library(Rmisc)
library(cowplot)

dir <- "~/Dropbox/dropbox_work/data/GWAS_m_f_fitness/"

###########################
## GWAS, male ####
###########################

## Male
massoc <- read.table(paste0(dir,"male_maf0.05.assoc"),head=T)
massoc$Chromosome <- as.factor(massoc$Chromosome)
levels(massoc$Chromosome) <- c("2L","2R","3L","3R","X")
#Q values
massoc$fdr <- p.adjust(massoc$Wald_P,"BH")

#Manhattan plot, without Bonferroni line
CMplot(massoc[,c("Predictor","Chromosome","Basepair","Wald_P")],type="p",plot.type="m",col=c("lightblue","dodgerblue4"),LOG10=TRUE,threshold=c(min(massoc$Wald_P[massoc$fdr>=0.3])),amplify = F,threshold.col=c("grey"),memo="",verbose=TRUE,width=14,height=6,cex=0.5,file="tiff",dpi=600,chr.labels = c("2L","2R","3L","3R","X"))
#Manhattan plot, with Bonferroni line
#CMplot(massoc[,c("Predictor","Chromosome","Basepair","Wald_P")],type="p",plot.type="m",col=c("lightblue","dodgerblue4"),LOG10=TRUE,threshold=c(min(massoc$Wald_P[massoc$fdr>=0.3]),0.05/nrow(massoc)),threshold.col=c("grey","black"),memo="",verbose=TRUE,width=14,height=6,cex=0.5,file="tiff",dpi=600,chr.labels = c("2L","2R","3L","3R","X"))
#CMplot(massoc[,c("Predictor","Chromosome","Basepair","Wald_P")],type="p",plot.type="q",col="dodgerblue4",LOG10=TRUE,threshold=NULL,memo="",dpi=300,verbose=TRUE,file="tiff",width=14,height=6)
names(massoc)[2] <- "SNP"

###########################
## GWAS, female ####
###########################

#Import LDAK mixed model association analysis
fassoc <- read.table(paste0(dir,"female_maf0.05.assoc"),head=T)
fassoc$Chromosome <- as.factor(fassoc$Chromosome)
levels(fassoc$Chromosome) <- c("2L","2R","3L","3R","X")
#Q values
fassoc$fdr <- p.adjust(fassoc$Wald_P,"BH")

#Manhattan plot, with Bonferroni line
CMplot(fassoc[,c("Predictor","Chromosome","Basepair","Wald_P")],type="p",plot.type="m",col=c("darksalmon","firebrick4"),LOG10=TRUE,amplify = F,threshold.col=c("grey"),memo="",verbose=TRUE,width=14,height=6,cex=0.5,file="tiff",dpi=600,chr.labels = c("2L","2R","3L","3R","X"))
#Manhattan plot, without Bonferroni line
#CMplot(fassoc[,c("Predictor","Chromosome","Basepair","Wald_P")],type="p",plot.type="m",col=c("darksalmon","firebrick4"),LOG10=TRUE,threshold=c(0.05/nrow(massoc)),threshold.col=c("black"),memo="",file="tiff",dpi=600,verbose=TRUE,width=14,height=6,cex=0.5,chr.labels = c("2L","2R","3L","3R","X"))
#CMplot(fassoc[,c("Predictor","Chromosome","Basepair","Wald_P")],type="p",plot.type="q",col="firebrick4",LOG10=TRUE,threshold=NULL,memo="",dpi=300,verbose=TRUE,file="tiff",width=14,height=6)
names(fassoc)[2] <- "SNP"

##Histogram of P/Q values
par(mfrow=c(2,1))
hist(fassoc$Wald_P,breaks=100,col=rgb(1,0,0,1/2),main="",xlab="p values")
hist(massoc$Wald_P,breaks=100,col=rgb(0,0,1,1/2),add=T)
hist(massoc$fdr,breaks=100,col=rgb(0,0,1,1/2),main="",xlab="q values")
hist(fassoc$fdr,breaks=50,col=rgb(1,0,0,1/2),add=T)

## Genomic inflation factors
#Males
chisq <- qchisq(1 - massoc$Wald_P, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
#Females
chisq <- qchisq(1 - fassoc$Wald_P, 1)
lambda <- median(chisq) / qchisq(0.5, 1)


###########################
## Permuted phenotypes ####
###########################

## Create 1000 'pheno.txt' files where individual labels have been shuffled 
#Other columns are otherwise unchanged
#Delta = -0.25
set.seed(123)
for (i in 1:1000){
  pheno.m <- read.table(paste0(dir,"reml_maf0.05/male_maf0.05_reml.indi.res"),head=T)
  pheno.f <- read.table(paste0(dir,"reml_maf0.05/female_maf0.05_reml.indi.res"),head=T)
  newcols <- sample(pheno.m[,1],replace=F,size=nrow(pheno.m))
  pheno.m$IID <- newcols
  pheno.m$FID <- newcols
  pheno <- cbind(pheno.m[,c("IID","FID","Residual")],pheno.f[,c("Residual")])
  names(pheno)[3:4] <- c("Residual_male_fitness","Residual_female_fitness")
  write.table(pheno,paste(paste0(dir,"1000_permuted_phenos/pheno"),i,".txt",sep = "_"),row.names=F,col.names=F,quote=F)
  print(i)
}

#Delta = -1
set.seed(123)
for (i in 1:1000){
  pheno.m <- read.table(paste0(dir,"reml_maf0.05/male_maf0.05_weight_minus_1_reml.indi.res"),head=T)
  pheno.f <- read.table(paste0(dir,"reml_maf0.05/female_maf0.05_weight_minus_1_reml.indi.res"),head=T)
  newcols <- sample(pheno.m[,1],replace=F,size=nrow(pheno.m))
  pheno.m$IID <- newcols
  pheno.m$FID <- newcols
  pheno <- cbind(pheno.m[,c("IID","FID","Residual")],pheno.f[,c("Residual")])
  names(pheno)[3:4] <- c("Residual_male_fitness","Residual_female_fitness")
  write.table(pheno,paste(paste0(dir,"1000_permuted_phenos_weight_minus_1/pheno"),i,".txt",sep = "_"),row.names=F,col.names=F,quote=F)
  print(i)
}

#Delta = 0
set.seed(123)
for (i in 1:1000){
  pheno.m <- read.table(paste0(dir,"reml_maf0.05/male_maf0.05_weight_0_reml.indi.res"),head=T)
  pheno.f <- read.table(paste0(dir,"reml_maf0.05/female_maf0.05_weight_0_reml.indi.res"),head=T)
  newcols <- sample(pheno.m[,1],replace=F,size=nrow(pheno.m))
  pheno.m$IID <- newcols
  pheno.m$FID <- newcols
  pheno <- cbind(pheno.m[,c("IID","FID","Residual")],pheno.f[,c("Residual")])
  names(pheno)[3:4] <- c("Residual_male_fitness","Residual_female_fitness")
  write.table(pheno,paste(paste0(dir,"1000_permuted_phenos_weight_0/pheno"),i,".txt",sep = "_"),row.names=F,col.names=F,quote=F)
  print(i)
}
