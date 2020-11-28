library(CMplot)
library(qqman)
library(Hmisc)
library(reshape2)
library(ggplot2)
library(data.table)
library(Rmisc)
library(cowplot)

dir <- "~/Dropbox/dropbox_work/data/GWAS_m_f_fitness/"

#Import GWAS data
massoc <- read.table(paste0(dir,"male_maf0.05.assoc"),head=T)
massoc$Chromosome <- as.factor(massoc$Chromosome)
levels(massoc$Chromosome) <- c("2L","2R","3L","3R","X")
names(massoc)[2] <- "SNP"

fassoc <- read.table(paste0(dir,"female_maf0.05.assoc"),head=T)
fassoc$Chromosome <- as.factor(fassoc$Chromosome)
levels(fassoc$Chromosome) <- c("2L","2R","3L","3R","X")
fassoc$fdr <- p.adjust(fassoc$Wald_P,"BH")
names(fassoc)[2] <- "SNP"


##########################
## SNP variant effects####
##########################

snp.effects <- read.table(paste0(dir,"VEP_all_severe.txt"))
snp.effects <- snp.effects$V4
fassoc$Variant_effect <- snp.effects
massoc$Variant_effect <- snp.effects

##########################
## LD pruning ####
##########################

#Randomised P-value column (for LD pruning)
set.seed(123)
fassoc$Wald_P_randomised <- sample(fassoc$Wald_P,size = nrow(fassoc),replace=F)
write.table(fassoc,paste0(dir,"female_maf0.05.assoc2"),quote=F,row.names=F)

## Randomised P-value column for coding sites only
write.table(subset(fassoc,Variant_effect %in% c("synonymous_variant","missense_variant")),paste0(dir,"female_maf0.05.assoc3"),quote=F,row.names=F)

## Import LD-pruned data

#See UNIX code for obtaining female_maf0.05_pruned_r0.4_10kb.clumped
fassoc.pruned <- read.table(paste0(dir,"female_maf0.05_pruned_r0.4_10kb.clumped"),head=T)
fassoc$LD_indep <- ifelse(fassoc$SNP %in% fassoc.pruned$SNP,1,0)
massoc$LD_indep <- ifelse(massoc$SNP %in% fassoc.pruned$SNP,1,0)
#See UNIX code for obtaining female_maf0.05_coding_pruned_r0.4_10kb.clumped
fassoc.pruned.coding <- read.table(paste0(dir,"female_maf0.05_coding_pruned_r0.4_10kb.clumped"),head=T)
fassoc$LD_indep_coding <- ifelse(fassoc$SNP %in% fassoc.pruned.coding$SNP,1,0)
massoc$LD_indep_coding <- ifelse(massoc$SNP %in% fassoc.pruned.coding$SNP,1,0)


##########################
## Candidates, genomic location ####
##########################

#Disproportionately X-linked or autoosmal?
m.clust <- read.table(paste0(dir,"clumped_highly_associated_r0.4_10kb.clumped"),head=T)
m.clust <- m.clust$SNP
obs.X <- 16
obs.A <- 15
exp.X <- sum(massoc$LD_indep==1 & massoc$Chromosome=='X')
exp.A <- sum(massoc$LD_indep==1 & massoc$Chromosome!='X')
chisq.test(cbind(c(obs.X,obs.A),c(exp.X,exp.A)))

## Candidate genes
gbat.f <- read.table(paste0(dir,"gbat_female_maf0.05/remls.1"),head=T)
gbat.f$fdr <- p.adjust(gbat.f$LRT_P_Raw,"BH")
gbat.f <- gbat.f[order(gbat.f$fdr),]
cand.f <- data.frame(subset(gbat.f,fdr<0.3)$Gene_Name)
names(cand.f) <- "Gene_Name"

gbat.m <- read.table(paste0(dir,"gbat_male_maf0.05/remls.1"),head=T)
gbat.m$fdr <- p.adjust(gbat.m$LRT_P_Raw,"BH")
gbat.m <- gbat.m[order(gbat.m$fdr),]
cand.m <- data.frame(subset(gbat.m,fdr<0.3)$Gene_Name)
names(cand.m) <- "Gene_Name"

#Candidate gene table
gbat.f2 <- gbat.f[gbat.f$fdr<0.3,c("Gene_Name","LRT_P_Raw","fdr")]
gbat.f2$Sex <- "Female"
gbat.f2 <- gbat.f2[1:70,]
gbat.m2 <- gbat.m[gbat.m$fdr<0.3,c("Gene_Name","LRT_P_Raw","fdr")]
gbat.m2$Sex <- "Male"
gbat.m2 <- gbat.m2[1:22,]
write.table(rbind(gbat.f2,gbat.m2),paste0(dir,"candidate_genes_table.txt"),quote=F,row.names=F,sep="\t")

#Are they enriched on the X?
ensgenes <- read.table(paste0(dir,"ensgenes.txt"))

#Females
cand.f.X <- nrow(subset(cand.f,Gene_Name %in% subset(ensgenes,V2==5)$V1))
cand.f.A <- nrow(subset(cand.f,Gene_Name %in% subset(ensgenes,V2!=5)$V1))
noncand.f.X <- nrow(subset(gbat.f,Gene_Name %in% subset(ensgenes,V2==5)$V1))
noncand.f.A <- nrow(subset(gbat.f,Gene_Name %in% subset(ensgenes,V2!=5)$V1))
chisq.test(cbind(c(cand.f.X,cand.f.A),c(noncand.f.X,noncand.f.A)))
chisq.test(cbind(c(cand.f.X,cand.f.A),c(noncand.f.X,noncand.f.A)))$expected
#OR
(9/61)/(2197/12815)

#Males
cand.m.X <- nrow(subset(cand.m,Gene_Name %in% subset(ensgenes,V2==5)$V1))
cand.m.A <- nrow(subset(cand.m,Gene_Name %in% subset(ensgenes,V2!=5)$V1))
noncand.m.X <- nrow(subset(gbat.m,Gene_Name %in% subset(ensgenes,V2==5)$V1))
noncand.m.A <- nrow(subset(gbat.m,Gene_Name %in% subset(ensgenes,V2!=5)$V1))
chisq.test(cbind(c(cand.m.X,cand.m.A),c(noncand.m.X,noncand.m.A)))
chisq.test(cbind(c(cand.m.X,cand.m.A),c(noncand.m.X,noncand.m.A)))$expected
#OR
(16/6)/(2197/12815)

##########################
## Candidates, variant effects ####
##########################

#SNP variant effects, male candidates
massoc$Candidate <- ifelse(massoc$SNP %in% m.clust,1,0)

pval <- vector()
obs_fun<- vector()
exp_fun<- vector()
or <- vector()
for (i in 1:16){
  obs.fun <- sum(massoc$Candidate==1 & massoc$Variant_effect==levels(massoc$Variant_effect)[i])
  obs.other <- sum(massoc$Candidate==1 & massoc$Variant_effect!=levels(massoc$Variant_effect)[i])
  exp.fun <- sum(massoc$LD_indep==1  & massoc$Variant_effect==levels(massoc$Variant_effect)[i])
  exp.other <- sum(massoc$LD_indep==1 & massoc$Variant_effect!=levels(massoc$Variant_effect)[i])
  pval[i] <- as.numeric(chisq.test(cbind(c(obs.fun,obs.other),c(exp.fun,exp.other)))[3])
  obs_fun[i] <- obs.fun
  exp_fun[i] <- exp.fun
  or[i] <- (obs.fun/obs.other)/(exp.fun/exp.other)
}
m.cand.vep <- data.frame(cbind(levels(massoc$Variant_effect),pval,obs_fun,exp_fun,or))
