#########################################################################################################
# Jackknife Conformal method (calibration)
#########################################################################################################
# Arguments:

# x: A matrix or a data frame whose rows contain the observed trajectories for a given subject.
# mesh: grid points at which the B-spline basis is evaluated.
# dimension.Bspline: number of B-spline basis functions.
# q: chosen dimension of the functional PCA
# alpha: miscoverage level of the prediction region. Thus, (1 - alpha) is the nominal coverage level.

#########################################################################################################

library(mgcv)

Jackknife.fun.calib <- function(x,
                          mesh, dimension.Bspline,
                          q, alpha)
{
# R functions to generate B-spline basis
  Bspline.manual <- function(mesh, dimension.Bspline){
    if (dimension.Bspline == 5)
      part <- 2
    else
      part <- dimension.Bspline/2
    aa <- ((-part):part)/part
    dimension.Bsplinenots <- (aa+1)*1.0/2
    base.estim.disc <- cSplineDes(mesh,dimension.Bsplinenots)
    base.ortonorm.disc <- qr.Q(qr(base.estim.disc))
    return(base.ortonorm.disc)
  }
  
# Generating B-spline basis functions
  base.estim.Bspline <- Bspline.manual(mesh, dimension.Bspline) 
  
# Function to compute standardized residuals as a weighted Euclidean distance
  distancescore <- function(y,z,w){
    out <- sqrt(sum(((y-z)^2) * w))
    return(out)
  }
  
  x<-data.matrix(x)
  
# Projecting subject data onto the B-spline basis
  y <- x %*% base.estim.Bspline
  n <- nrow(y)
  p <- ncol(y)
  
# Variance estimation for the standardization of residuals
var.explained <- sapply(1:p, function(k) {
  y.DsubLTS <- DsubLTS(x = y, q = k, alpha = alpha, scale=FALSE)
    mu.y.DsubLTS <- y.DsubLTS$mu
    bon <- y.DsubLTS$B
    yc <- y - matrix(mu.y.DsubLTS, nrow = n, ncol = p, byrow = TRUE)
    y.rulo.DsubLTS <- ((yc %*% bon) %*% t(bon)) + 
      matrix(mu.y.DsubLTS, nrow = n, ncol = p, byrow = TRUE)

    1 - (mean(rowSums((y- y.rulo.DsubLTS)^2)) /
           mean(rowSums((y - mu.y.DsubLTS)^2)))
  })
  
# Standardization weights
  pesos <- pmax(0, c(var.explained[1], diff(var.explained)))
  pesos <- pesos / sum(pesos);pesos
  

# Computing standardized Jackknife residuals on the calibration set
  resid.jackknife.pred <- sapply(1:n, function(i) {
    # Jackknife calibration using robust PCA (DsubLTS algorithm)
    jackknife.fit.DsubLTS <- DsubLTS(x = y[-i, ], q = q, alpha = alpha, scale=FALSE)
    yc <- y[i, ] - jackknife.fit.DsubLTS$mu
    bon <- jackknife.fit.DsubLTS$B
    jackknife.DsubLTS.i.prediction <- ((yc %*% bon) %*% t(bon)) + 
      jackknife.fit.DsubLTS$mu
    # Standardized Jackknife residuals on the calibration set
    distancescore(y[i, ], jackknife.DsubLTS.i.prediction, w = pesos)
  })
  
# Computing the $(1-\alpha)$-quantile of
# the empirical distribution of the calibration residuals.
  cuantil <- quantile(resid.jackknife.pred, probs= alpha)
  
# Calibration using robust PCA (DsubLTS algorithm)
  y.DsubLTS <- DsubLTS(x=y, q=q, alpha=alpha, scale=FALSE)
  mu.y.DsubLTS <- y.DsubLTS$mu
  bon <- y.DsubLTS$B
  
  return(list(
    quantile = cuantil,
    calib.proj.mu = mu.y.DsubLTS,
    calib.proj.B = bon,
    calib.proj.weights = pesos
  ))
}