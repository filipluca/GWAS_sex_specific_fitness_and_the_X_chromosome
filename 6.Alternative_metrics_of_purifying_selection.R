library(ggplot2)
library(cowplot)
library(reshape2)

dir <- "~/Dropbox/dropbox_work/data/GWAS_m_f_fitness/"

##############
## Pn, Ps ####
##############

## Coding sites 
massoc$Nonsyn_or_syn <- factor(ifelse(massoc$Variant_effect=="missense_variant","NS",ifelse(massoc$Variant_effect=="synonymous_variant","S",NA)))
levels(massoc$Nonsyn_or_syn) <- c(1,0)
massoc$Nonsyn_or_syn <- as.numeric(as.character(massoc$Nonsyn_or_syn))

#Statistical test (no interaction term)
mod <- glm(data=subset(massoc,LD_indep_coding==1),Nonsyn_or_syn~MAF+factor(Chromosome=="X"),family="binomial")
summary(mod)
#95% confidence intervals
exp(confint(mod))

#Statistical test (with interaction term)
mod2 <- glm(data=subset(massoc,LD_indep_coding==1),Nonsyn_or_syn~MAF*factor(Chromosome=="X"),family="binomial")
summary(mod2)
#95% confidence intervals
exp(confint(mod2))


#############
#Plots ##
#############

#Cut MAFs into bins (of equal frequency width, N=15)
massoc$MAF_bin <- factor(cut(massoc$MAF,breaks=seq(0.05,0.5,0.03)))
#Transform table for plotting
massoc.t <- data.frame(table(subset(massoc,LD_indep_coding==1)$Nonsyn_or_syn,subset(massoc,LD_indep_coding==1)$MAF_bin,subset(massoc,LD_indep_coding==1)$Chromosome=="X"))
names(massoc.t) <- c("Nonsyn_or_syn","MAF_bin","Is_X","Freq")
massoc.t$Nonsyn_or_syn <- as.factor(massoc.t$Nonsyn_or_syn)
levels(massoc.t$Nonsyn_or_syn) <- c("S","NS")
massoc.t2 <- dcast(massoc.t,Is_X+MAF_bin~Nonsyn_or_syn)
massoc.t2$NS_S_prop <- massoc.t2$NS/(massoc.t2$S+massoc.t2$NS)
levels(massoc.t2$MAF_bin) <- seq(0.05,0.5,0.03)+0.015
massoc.t2$MAF_bin <- as.numeric(as.character(massoc.t2$MAF_bin))

p2 <- ggplot(massoc.t2,aes(x=MAF_bin,y=NS_S_prop,col=Is_X))+
  geom_point(size=2)+
  theme_bw()+
  geom_smooth(method="lm")+
  ylim(c(0,0.4))+
  scale_color_manual(values=c("forestgreen","salmon"))+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),legend.position="none",strip.text=element_text(size=15),strip.background=element_rect(fill="white"))+
  xlab("MAF")+
  ylab("Pn / (Pn+Ps)")

p2

## Sup plots
par(mfrow=c(1,2))
hist(subset(massoc,LD_indep==1 & Chromosome!="X")$MAF,breaks=50,col=rgb(100,200,100,max = 255, alpha = 80),main="Autosomes",xlab="MAF",ylab="Count")
hist(subset(massoc,LD_indep==1 & Chromosome=="X")$MAF,breaks=50,col=rgb(200,100,50,max = 255, alpha = 80),main="X chromosome",xlab="MAF",ylab="Count")


