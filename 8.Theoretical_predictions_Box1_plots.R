##########################
## Box 1, Va plots
##########################

#r = sm/sf
F_over_M_fn_A <- function(h,r){
  (1/r)^2 }
F_over_M_fn_X <- function(h,r){
  2*h^2*(1/r)^2}
X_over_A_fn_F <- function(h,r){
  (3*h*(1+r)) / (2*(2*h+r)) }
X_over_A_fn_M <- function(h,r){
  (3*(1+r)) / (4*h*(2*h+(r))) }

va1a <- F_over_M_fn_A(rep(seq(0.02,1,by=0.02),2),c(rep(1,50),rep(2,50)))
va1b <- F_over_M_fn_X(rep(seq(0.02,1,by=0.02),2),c(rep(1,50),rep(2,50)))
va2a <- X_over_A_fn_F(rep(seq(0.02,1,by=0.02),2),c(rep(1,50),rep(2,50)))
va2b <- X_over_A_fn_M(rep(seq(0.02,1,by=0.02),2),c(rep(1,50),rep(2,50)))

par(mfrow=c(1,2))

#Autosomes, F/M
plot(rep(seq(0.02,1,by=0.02),2)[1:50],va1a[1:50],xlab="Dominance (h)",type="l",ylab=expression('V'["A,Female"] / 'V'["A,Male"]),ylim=c(0,3),col="forestgreen",lwd=2,main=expression('Female vs. male contribution to V'[A]))
points(0.5,F_over_M_fn_A(0.5,1),col="forestgreen",pch=19)
points(rep(seq(0.02,1,by=0.02),2)[51:100],va1a[51:100],type="l",col="forestgreen",lty=3,lwd=2)
points(0.5,F_over_M_fn_A(0.5,2),col="forestgreen")
abline(h=1,col="grey",lty=3)
#X, F/M
points(rep(seq(0.02,1,by=0.02),2)[1:50],va1b[1:50],type="l",ylim=c(0,3),col="darkorange",lwd=2)
points(0.5,F_over_M_fn_X(0.5,1),col="darkorange",pch=19)
points(rep(seq(0.02,1,by=0.02),2)[51:100],va1b[51:100],type="l",ylim=c(0,3),col="darkorange",lty=3,lwd=2)
points(0.5,F_over_M_fn_X(0.5,2),col="darkorange")

#Females, X/A
plot(rep(seq(0.02,1,by=0.02),2)[1:50],va2a[1:50],xlab="Dominance (h)",type="l",ylab=expression('V'["A,X"] / 'V'["A,Autosomes"]),ylim=c(0,3),col="red",lwd=2,main=expression('X vs. autosome contribution to V'[A]))
points(0.5,X_over_A_fn_F(0.5,1),col="red",pch=19)
points(rep(seq(0.02,1,by=0.02),2)[51:100],va2a[51:100],type="l",ylim=c(0,3),col="red",lty=3,lwd=2)
points(0.5,X_over_A_fn_F(0.5,2),col="red")
abline(h=1,col="grey",lty=3)
#Males, X/A
points(rep(seq(0.02,1,by=0.02),2)[1:50],va2b[1:50],xlab="h",type="l",ylim=c(0,3),col="blue",lwd=2)
points(0.5,X_over_A_fn_M(0.5,1),col="blue",pch=19)
points(rep(seq(0.02,1,by=0.02),2)[51:100],va2b[51:100],type="l",xlab="h",ylim=c(0,3),col="blue",lty=3,lwd=2)
points(0.5,X_over_A_fn_M(0.5,2),col="blue")


##########################
## Box 1, Regression of effect size on allele frequency plots
##########################

#r = sm/sf
beta_F_over_M_fn_A <- function(h,r){
  (1/r) }
beta_F_over_M_fn_X <- function(h,r){
  h*(1/r) }
beta_X_over_A_fn_F <- function(h,r){
  (2*(2*h+r)) / (3*h*(1+r)) }
beta_X_over_A_fn_M <- function(h,r){
  (2*(2*h+r)) / (3*h^2*(1+r)) }

reg1a <- beta_F_over_M_fn_A(rep(seq(0.02,1,by=0.02),2),c(rep(1,50),rep(2,50)))
reg1b <- beta_F_over_M_fn_X(rep(seq(0.02,1,by=0.02),2),c(rep(1,50),rep(2,50)))
reg2a <- beta_X_over_A_fn_F(rep(seq(0.02,1,by=0.02),2),c(rep(1,50),rep(2,50)))
reg2b <- beta_X_over_A_fn_M(rep(seq(0.02,1,by=0.02),2),c(rep(1,50),rep(2,50)))

par(mfrow=c(1,2))
#Autosomes, F/M
plot(rep(seq(0.02,1,by=0.02),2)[1:50],reg1a[1:50],xlab="Dominance (h)",type="l",ylab=expression(beta["Female"] / beta["Male"]),ylim=c(0,3),col="forestgreen",lwd=2,main=expression('Regression of effect size on allele frequency \n(Female vs. male)'))
points(0.5,beta_F_over_M_fn_A(0.5,1),col="forestgreen",pch=19)
points(rep(seq(0.02,1,by=0.02),2)[51:100],reg1a[51:100],type="l",col="forestgreen",lty=3,lwd=2)
points(0.5,beta_F_over_M_fn_A(0.5,2),col="forestgreen")
abline(h=1,col="grey",lty=3)
#X, F/M
points(rep(seq(0.02,1,by=0.02),2)[1:50],reg1b[1:50],type="l",ylim=c(0,3),col="darkorange",lwd=2)
points(0.5,beta_F_over_M_fn_X(0.5,1),col="darkorange",pch=19)
points(rep(seq(0.02,1,by=0.02),2)[51:100],reg1b[51:100],type="l",ylim=c(0,3),col="darkorange",lty=3,lwd=2)
points(0.5,beta_F_over_M_fn_X(0.5,2),col="darkorange")

#Females, X/A
plot(rep(seq(0.02,1,by=0.02),2)[1:50],reg2a[1:50],xlab="Dominance (h)",type="l",ylab=expression(beta["X"] / 'Beta'["Autosomes"]),ylim=c(0,3),col="red",lwd=2,main=expression('Regression of effect size on allele frequency \n(X vs. autosome)'))
points(0.5,beta_X_over_A_fn_F(0.5,1),col="red",pch=19)
points(rep(seq(0.02,1,by=0.02),2)[51:100],reg2a[51:100],type="l",ylim=c(0,3),col="red",lty=3,lwd=2)
points(0.5,beta_X_over_A_fn_F(0.5,2),col="red")
abline(h=1,col="grey",lty=3)
#Males, X/A
points(rep(seq(0.02,1,by=0.02),2)[1:50],reg2b[1:50],xlab="h",type="l",ylim=c(0,3),col="blue",lwd=2)
points(0.5,beta_X_over_A_fn_M(0.5,1),col="blue",pch=19)
points(rep(seq(0.02,1,by=0.02),2)[51:100],reg2b[51:100],type="l",xlab="h",ylim=c(0,3),col="blue",lty=3,lwd=2)
points(0.5,beta_X_over_A_fn_M(0.5,2),col="blue")
