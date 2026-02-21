library(mgcv)
library(EasyMMD)

# x: A list of subjects, where each element is a matrix whose rows contain the observed trajectories for a given subject.
# mesh: grid points at which the B-spline basis is evaluated.
# dimension.Bspline: number of B-spline basis functions.
# q: chosen dimension of the functional PCA
# indexTargetSubject: index of the element in x corresponding to the target subject.
# alpha: miscoverage level of the prediction region. Thus, (1 - alpha) is the nominal coverage level.
#####################################################################################

DistributionalShift.fun <- function(x,
                                    mesh, dimension.Bspline,
                                    q, 
                                    indexTargetSubject, alpha)
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

# Function to compute the $(1-\alpha)$-quantile of
# the weighted empirical distribution of the calibration residuals.
umbral.train <- function(x, w, p){
  sw <- sum(w)
  o <- order(x)
  xs <- x[o]
  ws <- w[o]
  cw <- cumsum(ws) / sw
  xs[ which(cw >= p)[1] ]
}

# Function to compute standardized residuals as a weighted Euclidean distance
distancescore <- function(y,z,w){
  out <- sqrt(sum(((y-z)^2) * w))
  return(out)
} 

list.scores.todos <- list() 
n.subjects <- length(x)

# Robust centering of projection data for each subject using the corresponding robust PCA center estimate
for (subject in 1:n.subjects) {
  x.obs <- x[[subject]] # Extract subject functional data

  # Projecting subject data onto the B-spline basis
  y <- x.obs %*% base.estim.Bspline
  n <- nrow(y)
  p <- ncol(y)
  
  # Centering projection data based on estimated center from DsubLTS algorithm
  y.DsubLTS <- DsubLTS(x=y, q=q, alpha=0.5, scale=FALSE)
  mu.y.DsubLTS <- y.DsubLTS$mu # Extracting estimated center
  ycc <- y - matrix(mu.y.DsubLTS, nrow = n, ncol = p, byrow = TRUE) # Robust centering
  list.scores.todos[[subject]] <- ycc
}
 
# Computing MMD matrix.
# Row index of the MMD matrix indicates target subject.
# Each row of the MMD matrix contains information of similarity of each of the other subjects
# with respect to that target subject.

dist_matrix <- matrix(0, nrow = n.subjects, ncol = n.subjects)

# Building the reference sequence
for (i.dmm in 1:n.subjects) {
  list_reference <- list.scores.todos
  tmp <- list_reference[[i.dmm]]
  list_reference[[i.dmm]] <- list_reference[[n.subjects]]
  list_reference[[n.subjects]] <- tmp
  
  matriz_A <- do.call(rbind, list_reference)
  
  # Building the swaped sequence
  for (j.dmm in 1:n.subjects) {
    
    list_j <- list_reference
    tmp <- list_j[[j.dmm]]
    list_j[[j.dmm]] <- list_j[[n.subjects]]
    list_j[[n.subjects]] <- tmp
    
    matriz_B <- do.call(rbind, list_j)
    
    xy <- rbind(t(matriz_A), t(matriz_B))
    
    # MMD distance with Gaussian kernel with a diagonal covariance matrix,
    # variance set to the squared median of the pairwise Euclidean distances 
    # between the feature vectors in the pooled data
    dists <- dist(xy)  
    median_dist <- median(as.numeric(dists))  
    var_kernel <- diag(rep(median_dist^2, ncol(t(matriz_A))))
  
    # Computing biased MMD distances
    dist_matrix[i.dmm, j.dmm] <- MMD(t(matriz_A), t(matriz_B), var = var_kernel, bias = TRUE)
  }
}

# Excluding subject indexTargetSubject for calibration
list.sin.indexTargetSubject <- list.scores.todos[-indexTargetSubject]

# Calibration set
y.sin.indexTargetSubject <- do.call(rbind, list.sin.indexTargetSubject)

################## Fitting the method on the calibration set ####################################

# The calibration set includes all curves from the other available individuals
n <- nrow(y.sin.indexTargetSubject)
p <- ncol(y.sin.indexTargetSubject)

# Variance estimation for the standardization of residuals
var.explained <- sapply(1:p, function(k) {
  y.DsubLTS <- DsubLTS(x = y.sin.indexTargetSubject, q = k, alpha = 0.5, scale=FALSE)
  mu.y.DsubLTS <- y.DsubLTS$mu
  yc <- y.sin.indexTargetSubject - matrix(mu.y.DsubLTS, nrow = n, ncol = p, byrow = TRUE)
  bon <- y.DsubLTS$B
  y.rulo.DsubLTS <- ((yc %*% bon) %*% t(bon)) + 
    matrix(mu.y.DsubLTS, nrow = n, ncol = p, byrow = TRUE)
  
  1 - (
    mean(rowSums((y.sin.indexTargetSubject - y.rulo.DsubLTS)^2)) /
      mean(rowSums((y.sin.indexTargetSubject - mu.y.DsubLTS)^2))
  )
})

# Standardization weights
pesos <- pmax(0, c(var.explained[1], diff(var.explained)))
pesos <- pesos / sum(pesos)

# Calibration using robust PCA (DsubLTS algorithm)
fit.DsubLTS <- DsubLTS(x = y.sin.indexTargetSubject, q = q, alpha = 0.5, scale=FALSE)
yc <- y.sin.indexTargetSubject - matrix(fit.DsubLTS$mu, nrow = n, ncol = p, byrow = TRUE)
bon <- fit.DsubLTS$B
DsubLTS.prediction <- ((yc %*% bon) %*% t(bon)) + 
  matrix(fit.DsubLTS$mu, nrow = n, ncol = p, byrow = TRUE)

# Residuals on the calibration set

# Creating names according to subject index
ids_originales <- seq_along(list.scores.todos) 
# Extracting index of the target individual
ids_filtrados <- ids_originales[-indexTargetSubject]  
# Repeating indexes of the other individuals according to the number of curves each of them has
paciente_id <- rep(ids_filtrados, times = sapply(list.sin.indexTargetSubject, nrow)) 

# Computing standardized residuals on the calibration set
resid.pred <- rep(0, n)
for(i.vir in 1:n){
  resid.pred[i.vir] <- distancescore(y.sin.indexTargetSubject[i.vir, ], 
                                     DsubLTS.prediction[i.vir,], 
                                     w = pesos)
  names(resid.pred)[i.vir] <- paste0("R", paciente_id[i.vir]) # assigning subject indices 
  # as names to the residuals
} 

# Computing weights for the other individuals relative to the target subject
# via the inverse-MMD rule

# The last swap leaves the sequence identical to the reference one (no-swap configuration), 
# as a result, the last column of the MMD distance matrix is equal to zero and is therefore removed.
D <- dist_matrix[ , -n.subjects] # MMD distance matrix
denominator <- colSums(1 / D)
W.subject.indexTargetSubject <- (1/D[indexTargetSubject,])/denominator

# Replicating weights according to the number of curves each of the other subjects have
# and assigning subject indices as names to the weights vector
nombres.R <- names(resid.pred) 
indices.R <- as.integer(sub("R", "", nombres.R)) 
W.rep <- W.subject.indexTargetSubject[indices.R]
indice.Rtarget <- which(nombres.R == paste0("R", n.subjects))
W.rep[indice.Rtarget] <- W.subject.indexTargetSubject[indexTargetSubject] 
names(W.rep) <- nombres.R

################## Predictions on the target subject’s projection data ####################################

# Computing predictions for the target subject’s projection data 
# using the robust PCA fit from the calibration set
y.nuevo <- list.scores.todos[[indexTargetSubject]]
n.y.nuevo <- dim(y.nuevo)[1]
yc <- y.nuevo - matrix(fit.DsubLTS$mu, nrow = n.y.nuevo, ncol = p, byrow = TRUE)
y.new.rulo.DsubLTS <- ((yc %*% bon) %*% t(bon)) + 
  matrix(fit.DsubLTS$mu, nrow = n.y.nuevo, ncol = p, byrow = TRUE)

# Computing standardized residuals on the projection data of the target subject
resid.y.new.pred <- sqrt(rowSums(((y.nuevo - y.new.rulo.DsubLTS)^2) *
                                   matrix(pesos, nrow = n.y.nuevo,
                                          ncol = p, byrow = TRUE)))

# Computing the $(1-\alpha)$-quantile of
# the weighted empirical distribution of the calibration residuals.
  cuantil <- umbral.train(unname(resid.pred), unname(W.rep), p=alpha) 

# Constructing the prediction region on the curves of the target subject
  index.good.new.curves <- which(resid.y.new.pred <= cuantil)
  
# Computing empirical coverage
cobertura.new.y <- length(index.good.new.curves)/n.y.nuevo

return(list(
  quantile = cuantil,
  predictionRegionIndices = index.good.new.curves,
  predictionRegionCurves = y.nuevo[index.good.new.curves, , drop = FALSE],
  coverage = cobertura.new.y
))
}