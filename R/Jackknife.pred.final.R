#########################################################################################################
# Jackknife Conformal method (prediction)
#########################################################################################################
# Arguments:

# xnew: A matrix or a data frame whose rows contain the new (future) trajectories for a given target subject.
# mesh: grid points at which the B-spline basis is evaluated.
# quantile: calibration quantile for a given miscoverage level of the prediction region.
# calib.proj.mu: jackknife estimated center obtained from the calibration projection data.
# calib.proj.B: jackknife estimated PCA directions obtained from the calibration projection data.
# calib.proj.weights: jackknife standardization weights obtained from the calibration projection data.

#########################################################################################################

library(mgcv)

Jackknife.fun.pred <- function(xnew,
                                mesh, dimension.Bspline,
                               quantile,
                               calib.proj.mu,
                               calib.proj.B,
                               calib.proj.weights)
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
  
  xnew <-data.matrix(xnew)
  
  # Projecting new data of the given target subject onto the B-spline basis
  y <- xnew %*% base.estim.Bspline
  n.y.nuevo <- nrow(y)
  p <- ncol(y)
  
  # Computing predictions on the projection data for the new trajectories of the given target subject 
  # using the robust PCA fit from the calibration set
  
  yc <- y - matrix(calib.proj.mu, nrow = n.y.nuevo, ncol = p, byrow = TRUE)
  y.new.rulo.DsubLTS <- ((yc %*% calib.proj.B) %*% t(calib.proj.B)) + 
    matrix(calib.proj.mu, nrow = n.y.nuevo, ncol = p, byrow = TRUE)
  
  resid.y.new.pred <- sqrt(rowSums(((y - y.new.rulo.DsubLTS)^2) * 
                                     matrix(calib.proj.weights, nrow = n.y.nuevo, 
                                            ncol = p, byrow = TRUE)))
  
  # Constructing the prediction region on the new curves of the target subject
  index.good.new.curves <- which(resid.y.new.pred <= quantile)
  
  # Computing empirical coverage
  coverage.new.y <- length(index.good.new.curves)/n.y.nuevo
  
  return(list(
    quantile = quantile,
    predictionRegionIndices = index.good.new.curves,
    predictionRegionCurves = xnew[index.good.new.curves, , drop = FALSE],
    coverage = coverage.new.y
  ))
}