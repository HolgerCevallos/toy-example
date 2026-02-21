#' @title Deterministic algorithm for the robust subspace LTS-estimator (DsubLTS)
#' @description  Computes the subspace LTS-estimator of Maronna (2005) using the fast and robust DsubLTS algorithm from Cevallos-Valdiviezo and Van Aelst (2019) which is an adaptation of the original algorithm of Maronna (2005). The DsubLTS algorithm directly estimates principal directions of the subspace without the need to calculate and decompose a covariance matrix, which makes it suitable for high-dimensional data. The algorithm starts by computing five robust deterministic values which are then improved iteratively by computing an iteratively reweighted least squares (IRLS) procedure in each iteration of the algorithm. Finally, the best solution is selected and iterated further until convergence. The IRLS procedure has been coded in C++ to make the algorithm even faster. 

#' @param x An (n x p) data matrix with observations in the rows and variables in the columns.
#' @param q dimension of the subspace. It must be less than or equal to the dimension of the data p.
#' @param N1 number of iterations to improve the scores A and the location mu, while keeping the initial loadings B fixed, for each deterministic value. Defaults to 3 as recommended by Maronna (2005).
#' @param N2 number of iterations to improve the scores A, the location mu and the initial loadings B for each deterministic value. Defaults to 2 as recommended by Maronna (2005).
#' @param N2bis number of iterations to further improve the estimates of A, B and mu corresponding to the best solution. Defaults to 10 as recommended by Maronna (2005).
#' @param Npc number of iterations in the IRLS procedure that updates A, B and mu in each iteration of the algorithm. Extensive experiments have shown that it suffices to use 3 iterations (default) to obtain stable results. See Details.
#' @param tol convergence criterion. The algorithm stops when successive changes of the objective function from one iteration to the next are less than or equal to \code{tol}. Default is 1e-06.
#' @param alpha trimming fraction used to compute the subspace LTS-estimator and to select clean subsets of the data in the computation of the deterministic starting values. Default is 0.5. See Details.
#' @param B.choice optional (p x q) matrix for the first q loadings. If the user provides the loadings b, the algorithm will only improve the corresponding estimates of the scores a and of the location mu. If a non-orthogonal matrix is provided, this function will first orthogonalize the provided loadings matrix before running the algorithm.

#' @return
#' \item{B}{(p x q) matrix containing the orthogonalized LTS-estimate for the first q loadings.}
#' \item{B.raw}{(p x q) matrix containing the raw LTS-estimate for the first q loadings.}
#' \item{mu}{p-dimensional vector containing the LTS-estimate for the location.}
#' \item{best.obj}{objective function value \eqn{\hat{\sigma}_{DsubLTS}(d(B.DsubLTS.raw,mu.DsubLTS))} corresponding to the final solution of the algorithm.}
#' \item{weights}{n-dimensional vector of (0-1) weights upon which the final subspace LTS-estimates were obtained. Observations with weight 1 are used to obtain the final robust estimates and could thus be identified as regular, while observations with weight 0 are not used and could thus be identified as potential outliers.}
#' \item{best.cand}{deterministic starting value strategy giving the final solution at the end of the algorithm. The algorithm uses five strategies to compute deterministic starting values: "rank", "tanh", "normalscore", "nontransformed" and "spatialsign". See Cevallos-Valdiviezo and Van Aelst (2019) for more details.}

#'@examples 
#' \dontrun{
#'# An example using the Octane data.
#'library(rrcov)
#'data(octane)
#'q <- 2 # subspace dimension
#'fit <- DsubLTS(x=octane, q=q, alpha=0.5)
#'fit$B.raw # raw B estimate
#'fit$B # orthogonalized B estimate
#' }

#' @references
#' \itemize{
#' \item{Cevallos-Valdiviezo, H., & Van Aelst, S. (2019). Fast computation of robust subspace estimators. Computational Statistics & Data Analysis, 134, 171-185.}
#' \item{Maronna, R. (2005). Principal components and orthogonal regression based on robust scales. Technometrics, 47(3), 264-273.}
#' }

#' @seealso S-L subspace estimator Cevallos-Valdiviezo and Van Aelst (2019) \code{\link{RsubLTS}}, Maronna(2005) \code{\link{SL}}.

#' @details The subspace LTS-estimator, as introduced by Maronna (2005), robustly estimates the best q-dimensional subspace by minimizing the LTS-scale of the Euclidean distances of the residuals. This function computes the subspace LTS-estimator using the fast and robust DsubLTS algorithm from Cevallos-Valdiviezo and Van Aelst (2019). Unlike the algorithm of Maronna (2005), the DsubLTS algorithm directly estimates the first q directions of the subspace without the need to calculate, store and decompose a covariance matrix. For that, DsubLTS uses an iterative reweighted least squares (IRLS) procedure to update solutions of the scores A, the loadings B and the location mu, in each iteration of the algorithm. The IRLS procedure makes use of the estimating equations of the estimator, which involve operations with q-dimensional vectors and matrices. Hence, the DsubLTS algorithm can be computed faster and is more suitable for high-dimensional data. In fact, it suffices to iterate the estimating equations only a few times (e.g. 3 times) to obtain close approximations to the first eigenvectors of the needed (weighted) covariance matrix. 
#' To reduce the computation time even further, the DsubLTS algorithm of Cevallos-Valdiviezo and Van Aelst (2019) uses five robust deterministic values to start the algorithm instead of using a large number of random orthogonal starting values. These five starting values can be computed fast as they involve an alternating least squares procedure requiring also q-dimensional vector and matrix operations. Hence, there is no need to compute a covariance matrix to obtain the initial values. Besides, these deterministic starting values uses simple transformations of the data and attempts to identify clean subsets of size \code{h=alpha*n} to obtain robust initial estimates. Thus, the algorithm starts close to a robust minimum of the objective function, which can generally result in quick convergence of the algorithm. 
#' This implementation can therefore be used to estimate the best q-dimensional subspace in the presence of outlying observations, specially in problems with high dimensions. For this deterministic algorithm, however, orthogonal equivariance is not guaranteed, although it has been shown empirically that in most cases there is not a serious loss of equivariance. Finally, this function outputs the raw final B estimate obtained at the end of the algorithm as well as a orthogonalized version of it. 

#' @import Rcpp
#' @import RcppArmadillo
#' @import robustbase
#' @import stats

#' @export

DsubLTS <- function(x, q, alpha=0.5, scale=TRUE, 
                    maxit1=3, maxit2=2, maxit2bis=10, maxitIRLS=3, tol=1e-6, 
                    B.choice=NULL,
                    na.action=match.fun(getOption("na.action"))){
  
  #************************************************************************
  #You need to install the "Rcpp", "RcppArmadillo", "robustbase" packages
  #require(Rcpp)
  #require(RcppArmadillo)
  #library(robustbase)
  
  # You need to load the following ".cpp" files:
  #sourceCpp('MyclassicalPcs.v2.cpp')
  #sourceCpp('MyPcs.cpp')
  #sourceCpp('MyPcs2.cpp')
  
  #************************************************************************
  x<-data.matrix(x) # NUEVO!
  
  if(anyNA(x)){
    if(missing(na.action)){
      naa <- getOption("na.action")
      if(is.null(naa)){
        na.action <- na.fail
        x <- na.action(x)
      }
      else if(is.character(naa)){
        na.action <- match.fun(na.action)
        x <- na.action(x)
      } 
    }else if(is.null(na.action)){ 
      na.action <- na.pass
      x <- na.action(x)
    }else{
      valid_actions <- list(na.fail, na.omit, na.exclude, na.pass)
      if(any(vapply(valid_actions, identical, logical(1), na.action))){
        x <- na.action(x)
      }else{
        stop("Invalid function specified in 'na.action' argument. It must be one of: na.fail, na.omit, na.exclude, or na.pass.")
      }
    }
  }
  
  n<-dim(x)[1]
  p<-dim(x)[2]
  
  if(alpha<0 | alpha>0.5){ # NUEVO!
    stop('Allowed values for alpha (trimming fraction) are between 0 and 0.5.')
  }
  
  if(q>p){ # NUEVO!
    q <- p-1 
    warning('The subspace dimension q is greater than p. The value of q has been set to p-1.') 
  } 
  
  h<-n-floor(n*alpha)
  keep_s <- Inf
  
  ##################################################
  # Computing scale for standardization
  ##################################################
  
  # Qnscales=apply(x,2,Qn) # NUEVO!
  
  # Compute scales of the variables if requested
  if (is.logical(scale)){
    if (!scale) scale.x <- rep(1, p) # vector("numeric", d) + 1
    if (scale){
      scale.x <- apply(x,2,Qn, na.rm = TRUE)
    }
  }else if(is.function(scale)){
    scale.x <- apply(x, 2, FUN = scale, na.rm = TRUE)
  }else if(is.null(scale)){
    scale.x <- rep(1, p)
  }else if(is.vector(scale)&&length(scale)==p&&is.numeric(scale)){
      scale.x <- scale
  }else{
    stop("'scale' must be logical, a function, numeric vector of length p, or NULL")
  }   
  
  # if (any(scale.x<1e-10)){
  #   small.scales.vars <- colnames(x)[which(scale.x < 1e-10)] 
  #   warning("The following variables have a scale too close to zero: ", 
  #           paste(small.scales.vars, collapse = ", "))
  # }
  
  if (any(scale.x<0.001)){
    small.scales.vars <- colnames(x)[scale.x < 0.001] 
    warning("The following variables have a scale too close to zero: ", 
            paste(small.scales.vars, collapse = ", "), 
            ". No scaling will be performed.")
    scale.x <- rep(1, p)
  }
  
  # Scale the variables
  x   <- x  / # NUEVO!
    matrix(scale.x,   nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  # x   <- sweep(x, 2, scale.x, "/")  # Note that scale.x is 1 when scale=F
  
  #####################################################
  
  # Give an initial (clean) h-subset
  determin.start <- function(x, q, scale=scale, #Qn.coorwise=Qnscales, # NUEVO!
                             strategy, alpha, Npc, tol){
    n <- dim(x)[1] # NUEVO!
    p <- dim(x)[2] # NUEVO!
    h <- n - floor(n*alpha)
    if(scale==FALSE | is.null(scale)){ # NUEVO!
      x.est <- (x - matrix(apply(x,2,median), nrow = n, ncol = p, byrow = TRUE)) /
        matrix(apply(x,2,Qn), nrow = n, ncol = p, byrow = TRUE) # NUEVO!
    }else{
      x.est <- x - rep(apply(x,2,median), rep.int(n, p))
    }
    # x.est <- scale( x , center= apply(x,2,median), 
    #                 scale = scaling) # NUEVO!
    if(strategy=="rank"){
      U <- apply(x, 2, rank)
    }else if(strategy=="tanh"){
      Y <- apply(x.est, 2, tanh)
      U <- t((1/apply(Y,2,Qn)) * t(Y))
    }else if(strategy=="normalscore"){
      U <- qnorm( (apply(x, 2, rank) - (1/3) ) / (n + (1/3) )  )
    }else if(strategy=="nontransformed"){
      U <- x.est
    }else if(strategy=="spatialsign"){
      Y <- (1/sqrt( apply(x.est^2, 1, sum))) * x.est
      U <- t((1/apply(Y,2,Qn)) * t(Y))
    }
    # Classical PCA on transformed matrix U
    mu <- apply(U,2,mean)
    b <- rbind(diag(q),matrix(0,nrow=p-q,ncol=q))
    update=findpcs0(x=U, B=b, mu=mu, Npc=Npc, tol=tol)
    # b=update$Bmat # NUEVO!
    # PC.directions <- qr.Q(qr(b))
    PC.directions <- update$Borthomat # NUEVO!
    
    Z <- crossprod(t(x.est), PC.directions)
    # Z <- scale( U , center= apply(U,2,median), scale= apply(U,2,Qn))
    obswise.sqnorm <- apply(Z^2, 1, sum)
    U <- x[obswise.sqnorm <= sort(obswise.sqnorm)[h],]
    return(U)
  }
  
  cand <- c("rank", "tanh", "normalscore", "nontransformed", "spatialsign") # NUEVO!
  
  if(is.null(B.choice)){ # NUEVO!
  
    for(i in cand){
      # Deterministic starting value for B
      U <- determin.start(x=x, q=q, scale=scale, #Qn.coorwise=Qnscales, # NUEVO!
                          strategy=i, alpha=alpha,
                          Npc=maxitIRLS, tol=tol)
      mu <- apply(U,2,mean)
      b <- rbind(diag(q),matrix(0,nrow=p-q,ncol=q))
      update=findpcs0(x=U, B=b, mu=mu, Npc=maxitIRLS, tol=tol)
      b=update$Bmat
      mu=update$mu
      
      update=findpcs2(x=x, B=b, mu=mu, h=h, s=0, N1=maxit1, N2=maxit2, Npc=maxitIRLS, tol=tol,
                      fixed=FALSE) 
      b=update$Bmat
      mu=update$mu
      s=update$scale

      # Storage of the b and mu with the lowest value of s
      if(s < keep_s){
        opt_cand <- i
        keep_s <- s
        keep_mu <- mu
        keep_b <- b
      }
    } # End of loop over the five deterministic starting values
    
    mu <- keep_mu
    b <- keep_b
    s <- keep_s
    
    # Improvement of the best deterministic starting value
    update=findpcs(x=x, B=b, mu=mu, h=h, s=s, N1=0, N2=maxit2bis, Npc=maxitIRLS, tol=tol,
                   fixed=FALSE)  
    b=update$Borthomat # NUEVO!
    b.raw=update$Bmat # NUEVO!
    mu=update$mu
    s=update$scale
    best <- update$index1 # NUEVO!
    # ws<-rep(0,n) # NUEVO!
    # ws[update$index2] <-1 # NUEVO!
    scores=update$Amat
    
  }else{
    if(!is.matrix(B.choice) | nrow(B.choice)!= p | 
       ncol(B.choice) != q | !is.numeric(B.choice)){
      stop(paste('B.choice must be a numeric matrix of size', p,'times', q,'.') ) # NUEVO!
    }else{ # NUEVO!
      b <- qr.Q(qr(B.choice)) # NUEVO!
      
      for(i in cand){
        
        # Deterministic starting value for B
        U <- determin.start(x=x, q=q, scale=scale, #Qn.coorwise=Qnscales, # NUEVO!
                            strategy=i, alpha=alpha,
                            Npc=maxitIRLS, tol=tol)
        mu <- apply(U,2,mean)
        
        update=findpcs2(x=x, B=b, mu=mu, h=h, s=0, N1=maxit1, N2=maxit2, Npc=maxitIRLS, tol=tol,
                        fixed=TRUE) 
        b=update$Bmat
        mu=update$mu
        s=update$scale

        # Storage of the b and mu with the lowest value of s
        if(s < keep_s){
          opt_cand <- i
          keep_s <- s
          keep_mu <- mu
          keep_b <- b
        }
      } # End of loop over the five deterministic starting values
      
      mu <- keep_mu
      b <- keep_b
      s <- keep_s
      
      # Improvement of the best deterministic starting value
      update=findpcs(x=x, B=b, mu=mu, h=h, s=s, N1=0, N2=maxit2bis, Npc=maxitIRLS, tol=tol,
                     fixed=TRUE)  
      b=update$Borthomat # NUEVO!
      b.raw=update$Bmat # NUEVO!
      mu=update$mu
      s=update$scale
      best <- update$index1 # NUEVO!
      # ws<-rep(0,n) # NUEVO!
      # ws[update$index2] <-1 # NUEVO!
      scores=update$Amat
    }
  }
  
  # Results
  res<-list(B.raw=b.raw, B=b, mu=as.vector(mu), 
            q=q, crit=s,  
            X=x, n=n, scores=scores, scale.x=scale.x, 
            alpha=alpha, h=h, best=as.vector(best), 
            iBest=opt_cand,
            na.action=na.action) # NUEVO!
  class(res) <- "DsubLTS"
  return(res)
}