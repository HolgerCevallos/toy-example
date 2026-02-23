require(Rcpp)
require(RcppArmadillo)
require(robustbase)
library(mgcv)
#install.packages("EasyMMD")
# install.packages("devtools")
# library(devtools)
# devtools::install_github("AnthonyEbert/EasyMMD")
library(EasyMMD)

# You need to load the following ".cpp" files to implement the robust DsubLTS algorithm:
sourceCpp('MyPcs.cpp')
sourceCpp('MyPcs2.cpp')
sourceCpp('MyclassicalPcs.v2.cpp')

# R code to implement the DsubLTS algorithm (robust PCA)
source("DsubLTS.R")

# R code to implement the Nonasymptotic Distributional Shift conformal method
source("NonasymptoticDistributionalShift.final.R")
source("Jackknife.calib.final.R")
source("Jackknife.pred.final.R")

###########################################################
# Design 1
###########################################################

lambda.k1.values <- function(k){
  lambda.k1 <- 0.5^(k-1)
  return(lambda.k1)
}

lambda.l2.values <- function(l){
  lambda.l2 <- 0.5^(l-1)
  return(lambda.l2)
}

xik.values <- function(lambda.k1){
  xik <- rnorm(1,0,lambda.k1)
  return(xik)
}

zetaijl.values <- function(lambda.l2){
  zetaijl <- rnorm(1,0,lambda.l2)
  return(zetaijl)
}

phik1.values <- function(t){
  phik1 <- list(sqrt(2)*sin(2*pi*t),
                sqrt(2)*cos(2*pi*t),
                sqrt(2)*sin(4*pi*t),
                sqrt(2)*cos(4*pi*t))
  return(phik1)
}

phil2.values <- function(t){
  phil2 <- list(sqrt(2)*sin(6*pi*t),
                sqrt(2)*cos(6*pi*t),
                sqrt(2)*sin(8*pi*t),
                sqrt(2)*cos(8*pi*t))
  return(phil2)
}

# epsilonij <- function(t, sigma){
#   set.seed(t)
#   rnorm(1, 0, sigma)
# }

##############################################################
# Data generation of calibration set
##############################################################

# Calibration set:
# 5 subjects, 50 curves are generated per subject at 21 time points

Npatients <- 5 #número de pacientes
Yij.todos <- vector("list", Npatients) 

for (s in 1:Npatients){
  
  sigma <- 1
  Nt <- 21
  curves.per.patient <- 50
  fixed.subject.effect <- rep(0,Nt)
  visit.effect <- matrix(0, nrow=curves.per.patient, ncol=Nt)
  Yij <- matrix(0, nrow=curves.per.patient, ncol=Nt)
  tm <- (0:20)/100
  
  for(t in 1:Nt){
    for(k in 1:4){
      lambda.k1 <- lambda.k1.values(k)
      set.seed(123+40+s)
      fixed.xk <- xik.values(lambda.k1)
      fixed.subject.effect[t] <- fixed.subject.effect[t] + 
        fixed.xk*phik1.values(tm[t])[[k]]
    }
  }
  
  for(j in 1:curves.per.patient){
    for(t in 1:Nt){
      for(l in 1:4){
        lambda.l2 <- lambda.l2.values(l)
        set.seed(123+s+j)
        zetaijl <- zetaijl.values(lambda.l2)
        visit.effect[j,t] <- visit.effect[j,t] + 
          zetaijl*phil2.values(tm[t])[[l]]
      }
    }
  }
  
  for(j in 1:curves.per.patient){
    for(t in 1:Nt){
      Yij[j,t] <- fixed.subject.effect[t] + visit.effect[j,t] +
        rnorm(1,0,sigma)
    }
  }
  Yij.todos[[s]] <- Yij
}

###################################################################################################
# Application of the Nonasymptotic Distributional Shift method to the toy data
###################################################################################################

# Target subject: subject 1
fit.pooled.method <- DistributionalShift.fun(x=Yij.todos,
                                    mesh=tm, dimension.Bspline=20,
                                    q=2, 
                                    indexTargetSubject=1, alpha=0.5)
fit.pooled.method$quantile
fit.pooled.method$predictionRegionIndices
fit.pooled.method$predictionRegionCurves
fit.pooled.method$coverage

###################################################################################################
# Application of the Jackknife method to the toy data
###################################################################################################

# Target subject: subject 1
fit.jackknife.method <- Jackknife.fun.calib(x=Yij.todos[[1]],
                                            mesh=tm, dimension.Bspline=20,
                                            q=2, alpha=0.5)

# Generation of 50 new curves for subject 1
s <- 1
sigma <- 1
Nt <- 21
curves.per.patient <- 50
fixed.subject.effect <- rep(0,Nt)
visit.effect <- matrix(0, nrow=curves.per.patient, ncol=Nt)
Yij <- matrix(0, nrow=curves.per.patient, ncol=Nt)
tm <- (0:20)/100

for(t in 1:Nt){
  for(k in 1:4){
    lambda.k1 <- lambda.k1.values(k)
    set.seed(1234+40+s)
    fixed.xk <- xik.values(lambda.k1)
    fixed.subject.effect[t] <- fixed.subject.effect[t] + 
      fixed.xk*phik1.values(tm[t])[[k]]
  }
}

for(j in 1:curves.per.patient){
  for(t in 1:Nt){
    for(l in 1:4){
      lambda.l2 <- lambda.l2.values(l)
      set.seed(1234+s+j)
      zetaijl <- zetaijl.values(lambda.l2)
      visit.effect[j,t] <- visit.effect[j,t] + 
        zetaijl*phil2.values(tm[t])[[l]]
    }
  }
}

for(j in 1:curves.per.patient){
  for(t in 1:Nt){
    Yij[j,t] <- fixed.subject.effect[t] + visit.effect[j,t] +
      rnorm(1,0,sigma)
  }
}
Yij.new.subject.1 <- Yij

pred.jackknife.method <- Jackknife.fun.pred(xnew=Yij.new.subject.1,
                  mesh=tm, dimension.Bspline=20,
                  quantile=fit.jackknife.method$quantile,
                  calib.proj.mu=fit.jackknife.method$calib.proj.mu,
                  calib.proj.B=fit.jackknife.method$calib.proj.B,
                  calib.proj.weights=fit.jackknife.method$calib.proj.weights)

pred.jackknife.method$quantile
pred.jackknife.method$predictionRegionIndices
pred.jackknife.method$predictionRegionCurves
pred.jackknife.method$coverage