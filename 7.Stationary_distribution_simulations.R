library(tidyverse)
library(devtools)
library(ggplot2)
library(GGally)
library(doParallel)
registerDoParallel(cores=4)

lower_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    scale_fill_viridis_b(alpha=0.7,direction = -1)+    
    stat_density2d(aes(fill=..density..), geom="tile", contour = FALSE)+
    geom_point(size=0.3)+
    stat_smooth(method="lm",col="black",se=F,size=0.7)+
    theme_bw()
  p
}
diag_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_histogram(bins=30,fill="grey", col="black") +
    theme_bw()
  p
}

##########################################
## Run stationary distribution simulations (test run, on 40 total simulations)
##########################################

output <- foreach(icount(40), .combine=rbind) %dopar% {
  
  #Mutation rate
  u = 10^-8
  #The DFE for non-neutral sites is gamma distributed, with shape and scale parameters k and theta (respectively)
  k = runif(1,min=0.25,max=0.4)
  s.avg = runif(1,min=10^-5,max=3.5*10^-3)
  theta = s.avg/k
  #Correlation between s.m and s.f
  correlation = runif(1)
  #Dominance 
  h = runif(1)
  #Effective population sizes
  Ne.A = runif(1,min=10^2,max=10^4)
  f.Ne.X = runif(1,min=0.5,max=1)
  Ne.X = Ne.A*f.Ne.X
  #Number of chromosomes sampled and sequenced
  n = 202
  #Number of sites per chromosome (synonymous and nonsynonymous) for autosomes
  L.synon.A = 10^7
  L.nonsynon.A = L.synon.A*2.5
  #Number of sites per chromosome (synonymous and nonsynonymous) for the X
  L.synon.X = L.synon.A*0.177
  L.nonsynon.X = L.nonsynon.A*0.177
  
  ###################
  ## Synon sites
  ###################
  
  #Simulated population-wide allele frequencies for each site
  freq.synon.A = rbeta(L.synon.A, 2*Ne.A*u, 2*Ne.A*u)
  freq.synon.X = rbeta(L.synon.X, 2*Ne.X*u, 2*Ne.X*u)
  #Checks
  #mean(freq.synon.A*(1-freq.synon.A)*2)
  #(2*Ne.A*u)/(1+4*Ne.A*u)
  #mean(freq.synon.X*(1-freq.synon.X)*2)
  #(2*Ne.X*u)/(1+4*Ne.X*u)
  
  ## Autosomes
  #vector to record the number of minor alleles at each site in the sample of sequenced chromosomes
  alleles.synon.A = rep(0, L.synon.A)
  #Simulation of minor allele counts in the sample of sequenced chromosomes
  alleles.synon.A <- rbinom(L.synon.A,n,freq.synon.A)
  #Folded
  alleles.synon.A <- ifelse(alleles.synon.A<(n/2),alleles.synon.A,n-alleles.synon.A)
  #segregating sites only, MAC>0.05
  alleles.synon.A = alleles.synon.A[alleles.synon.A>10]
  
  ## X chromosome
  alleles.synon.X = rep(0, L.synon.X)
  #Simulation of minor allele counts in the sample of sequenced chromosomes
  alleles.synon.X <- rbinom(L.synon.X,n,freq.synon.X)
  #Folded
  alleles.synon.X <- ifelse(alleles.synon.X<(n/2),alleles.synon.X,n-alleles.synon.X)
  #segregating sites only, MAC>0.05
  alleles.synon.X = alleles.synon.X[alleles.synon.X>10]
  
  ###################
  ## Nonsynon sites
  ###################
  
  #Simulation of selection coefficients for nonsynon sites
  sel.f = rep(0, L.nonsynon.A)
  sel.m <- rgamma(L.nonsynon.A, k, 1/theta)
  for(i in 1:(L.nonsynon.A)){
    if(abs(correlation) < 1){
      s.f <- rgamma(1, (k + rpois(1, sel.m[i]*correlation/(theta*(1-correlation)))), 1/(theta*(1-correlation)))
    } else{s.f <- correlation*sel.m[i]}
    sel.f[i] = s.f
  }
  
  ## Autosomes
  #Vector to hold population-wide allele frequencies for each site
  freq.nonsynon.A = rep(-1, L.nonsynon.A)
  x = runif(L.nonsynon.A)
  #Simulated population-wide allele frequencies for each autosomal site
  for(i in 1:L.nonsynon.A){
    #simulated allele frequency based on Smith and Connallon's (2017) rejection sampler
    while(freq.nonsynon.A[i]<0){
      y = rbeta(1, 2*Ne.A*u, 2*Ne.A*u)
      if(x[i] < exp(-(1/2)*Ne.A*(sel.m[i]+sel.f[i])*y*(2*h + y*(1 - 2*h)))){freq.nonsynon.A[i] = y} else{freq.nonsynon.A[i] = -1}
    }
  }
  #Checks
  #u/(s.avg*h)
  #Simulation of minor allele counts in the sample of sequenced chromosomes
  alleles.nonsynon.A <- rbinom(L.nonsynon.A,n,freq.nonsynon.A)
  #Folded
  alleles.nonsynon.A <- ifelse(alleles.nonsynon.A<(n/2),alleles.nonsynon.A,n-alleles.nonsynon.A)
  #MAC>0.05
  alleles.nonsynon.A = alleles.nonsynon.A[alleles.nonsynon.A>10]
  
  ## X chromosome
  #Vector to hold population-wide allele frequencies for each site
  freq.nonsynon.X = rep(-1, L.nonsynon.X)
  x = runif(L.nonsynon.X)
  #Simulated population-wide allele frequencies for each X-linked site
  for(i in 1:L.nonsynon.X){
    #simulated allele frequency based on Smith and Connallon's (2017) rejection sampler
    while(freq.nonsynon.X[i]<0){
      y = rbeta(1, 2*Ne.X*u, 2*Ne.X*u)
      if(x[i] < exp(-2*Ne.X*y*(sel.f[i]*(2*h + y*(1 - 2*h)) + sel.m[i])/3)){freq.nonsynon.X[i] = y} else{freq.nonsynon.X[i] = -1}
    }
  }
  #Checks
  #3*u/(2*s.avg*h+s.avg)
  #simulation of minor allele counts in the sample of sequenced chromosomes
  alleles.nonsynon.X <- rbinom(L.nonsynon.X,n,freq.nonsynon.X)
  #Folded
  alleles.nonsynon.X <- ifelse(alleles.nonsynon.X<(n/2),alleles.nonsynon.X,n-alleles.nonsynon.X)
  #MAC>0.05
  alleles.nonsynon.X = alleles.nonsynon.X[alleles.nonsynon.X>10]
  
  ## Dataframe with all necessary info to test conditions
  #Synon sites
  dd.neu1 <- data.frame(alleles.synon.A/n)
  names(dd.neu1) <- "MAF"
  dd.neu1$Chromosome <- "A"
  dd.neu1$Type <- 0
  dd.neu2 <- data.frame(alleles.synon.X/n)
  names(dd.neu2) <- "MAF"
  dd.neu2$Chromosome <- "X"
  dd.neu2$Type <- 0
  #Nonsynon sites
  dd.sel1 <- data.frame(alleles.nonsynon.A/n)
  names(dd.sel1) <- "MAF"
  dd.sel1$Chromosome <- "A"
  dd.sel1$Type <- 1
  dd.sel2 <- data.frame(alleles.nonsynon.X/n)
  names(dd.sel2) <- "MAF"
  dd.sel2$Chromosome <- "X"
  dd.sel2$Type <- 1
  #Both
  dd <- rbind(dd.neu1,dd.neu2,dd.sel1,dd.sel2)
  dd$Chromosome <- as.factor(dd$Chromosome)
  
  #Logistic regression
  X_effect <- exp(summary(glm(data=dd,Type~MAF*factor(Chromosome=="X"),family="binomial"))$coef[3,1])
  MAF_effect <- exp(summary(glm(data=dd,Type~MAF*factor(Chromosome=="X"),family="binomial"))$coef[2,1])
  Interaction_effect <- exp(summary(glm(data=dd,Type~MAF*factor(Chromosome=="X"),family="binomial"))$coef[4,1])
  
  #Outputs
  c(Ne.A,Ne.X,k,s.avg,correlation,h,X_effect,MAF_effect,Interaction_effect)
}

print(output)

##########################################
## Stationary distribution plots
##########################################

dir <- "~/Dropbox/dropbox_work/data/GWAS_m_f_fitness/"

#Clean output from cluster

file.names <- c("slurm-16616381.out","slurm-16616382.out","slurm-16616383.out","slurm-16616384.out","slurm-16616386.out","slurm-16616387.out","slurm-16616388.out","slurm-16616389.out","slurm-16616390.out","slurm-16616391.out")

dd2 <- list()
for (i in 1:length(file.names)){
  dd <- read.table(paste0(dir,"simulations/November/",file.names[i]),fill=T,skip=4)
  dd.tmp <- dd[1:2000,]
  dd.tmp2 <- dd[2002:4001,]
  dd <- cbind(dd.tmp,dd.tmp2[,2:5])
  for (j in 2:9){
    dd[,j] <- as.numeric(as.character(dd[,j]))
  }
  names(dd) <- c("sim","Ne.A", "Ne.X", "k", "s.avg", "cor.mf","h","X_coefficient","MAF_coefficient","Interaction_coefficient")
  dd2[[i]] <- dd
}

dd3 <- do.call(rbind,dd2)

dd4a <- subset(dd3,(X_coefficient>0.6273844 & X_coefficient<=0.9512530) & (MAF_coefficient>0.3441402 & MAF_coefficient<=0.6159652) & (Interaction_coefficient>0.5349367 & Interaction_coefficient<=2.5276752)  )[1:1000,]
dd4a$f.Ne.A <- dd4a$Ne.X/dd4a$Ne.A


#Median and quantiles 
quantile(dd4a$Ne.A,probs = c(0.025,0.5,0.975))
hist(dd4a$Ne.A)
quantile(dd4a$f.Ne.A,probs = c(0.025,0.5,0.975))
hist(dd4a$f.Ne.A)
quantile(dd4a$k,probs = c(0.025,0.5,0.975))
hist(dd4a$k)
quantile(dd4a$s.avg,probs = c(0.025,0.5,0.975))
hist(dd4a$s.avg)
quantile(dd4a$cor.mf,probs = c(0.025,0.5,0.975))
hist(dd4a$cor.mf)
quantile(dd4a$h,probs = c(0.025,0.5,0.975))
hist(dd4a$h)
quantile(dd4a$Ne.A*dd4a$s.avg*0.5,probs = c(0.025,0.5,0.975))
hist(dd4a$Ne.A*dd4a$s.avg*0.5)
quantile(dd4a$Ne.X*dd4a$s.avg*0.5,probs = c(0.025,0.5,0.975))
hist(dd4a$Ne.X*dd4a$s.avg*0.5)

#Correlation matrices
cor(dd4a[c("Ne.A","f.Ne.A","k","s.avg","cor.mf","h")],method = "spearman")
cor.test(dd4a$cor.mf,dd4a$h,method = "spearman")

## Pairs plot, Main 
ggpairs(dd4a[c("h","f.Ne.A")], lower=list(continuous=lower_fn),diag=list(continuous=diag_fn),upper=list(continuous="blank"),switch = 'y')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7,size=10),axis.text.y = element_text(size=10),panel.grid.major = element_blank(),strip.background =element_rect(fill="grey95"),strip.text.x =element_text(size=20),strip.text.y =element_text(size=10),panel.spacing = unit(1,"lines"))
## Pairs plot, Sup 
ggpairs(dd4a[c(2,11,4:7)], lower=list(continuous=lower_fn),diag=list(continuous=diag_fn),upper=list(continuous="cor",corMethod="spearman"),switch = 'y')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7,size=10),axis.text.y = element_text(size=10),panel.grid.major = element_blank(),strip.background =element_rect(fill="grey95"),strip.text.x =element_text(size=20),strip.text.y =element_text(size=10),panel.spacing = unit(1,"lines"))


