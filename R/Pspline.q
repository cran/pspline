smooth.Pspline <- function(x, y, w=rep(1, length(x)), norder=2, 
                           df=norder+2, spar=0, method=1) 
{

  #  Computes order NORDER polynomial smoothing spline:  the spline
  #    is piecewise of degree 2*NORDER - 1 and the norm of the
  #    derivative of order is penalized

  #  calls Fortran function spl

  #  Arguments:
  #  X       ...  argument values
  #  Y       ...  N by NVAR matrix of function values to be smoothed
  #  W       ...  weights (default all one's)
  #  NORDER  ...  order of smoothing spline (default 2)
  #  SPAR    ...  penalty parameter (default 0)
  #  DF      ...  effective degrees of freedom (trace(hatmatrix))
  #  METHOD  ...  smoothing method:  1  ...  fixed value of SPAR (default)
  #                                  2  ...  fixed value of DF
  #                                  3  ...  SPAR optimizes GCV criterion
  #                                  4  ...  SPAR optimizes  CV criterion

  #  Returns:  An object of class "smooth.Pspline" containing:
  #  NORDER  ...  order of smoothing spline
  #  X       ...  argument values
  #  YSMTH   ...  N by NVAR matrix of values of the smoothed functions
  #  LEV     ...  array of N leverage values
  #  GCV     ...  generalized cross-validation coefficient
  #  CV      ...  cross-validation coefficient
  #  DF      ...  final effective degrees of freedom
  #  SPAR    ...  final smoothing parameter value
  #  MY.CALL ...  calling statement

  my.call <- match.call()

  n <- length(x)
  if (is.matrix(y)) nvar <- ncol(y) else {
    nvar <- 1
    y <- as.matrix(y)
  }
  if (length(w) == 1) w <- rep(w,n)
  if (nrow(y) != n | length(w) != n) stop("Argument arrays of wrong length")
  if (method != 1 & method != 2 & method != 3 & method != 4) stop(
         "Wrong value for METHOD")
  if (norder <= 1 | norder >= 19) stop("Wrong value for NORDER")

  yhat     <- matrix(0,n,nvar)
  nworksiz <- (n-norder)*(4*norder + 3) + n 
  work     <- rep(0,nworksiz)
  lev      <- rep(0,n)
  gcv      <- 0
  cv       <- 0
  dfmax    <- n
  ier      <- 0
  irerun   <- 0

  result <- .Fortran("pspline", 
              as.integer(n),      as.integer(nvar),  as.integer(norder), 
              as.double(x),       as.double(w), 
              as.double(y),       as.double(yhat),   as.double(lev),
              as.double(gcv),     as.double(cv),     as.double(df), 
              as.double(spar),    as.double(dfmax), 
              as.double(work),    as.integer(method), 
              as.integer(irerun), as.integer(ier) )

  ier <- result[[17]]
  if (ier == 1) stop ("N < 2*NORDER + 1")
  if (ier == 2) stop ("NORDER < 2 or NORDER > 10")
  if (ier == 3) stop ("NVAR < 1")
  if (ier == 4) stop ("SPAR < 0")
  if (ier == 5) stop ("X not strictly increasing")
  if (ier == 6) stop ("W contains nonpositive values")
  if (ier < 0 ) stop ("Singularity error in solving equations")

  ysmth  <- matrix(result[[7]],n,nvar)
  lev    <- result[[8]]
  gcv    <- result[[9]]
  cv     <- result[[10]]
  df     <- result[[11]]
  spar   <- result[[12]]
  object <- list(norder = norder, x = x, ysmth = ysmth,  lev = lev, 
	         gcv = gcv, cv = cv,  df = df, 
		 spar = spar, call = my.call)
  class(object) <- "smooth.Pspline"
  object

}

     
  
  
predict.smooth.Pspline <- function(splobj, xarg, nderiv = 0) {
	
  if(missing(xarg)) return(splobj[c("x", "ysmth")])

  x      <- splobj$x
  ysmth  <- splobj$ysmth
  norder <- 2*splobj$norder
  n    <- length(x)
  nvar <- ncol(ysmth)
  narg <- length(xarg)

  if (nderiv < 0 | nderiv >= norder) stop("Violation of NDERIV >= NORDER.")
 
  dy    <- matrix(0,narg,nvar)
  work  <- rep(0,(2*norder+2)*n + norder)
  ier   <- 0

  result <- .Fortran("splifit", 
                     as.integer(n),    as.integer(narg),   
                     as.integer(nvar), as.integer(norder), as.integer(nderiv),
                     as.double(x),     as.double(ysmth),    
                     as.double(xarg),  as.double(dy),
                     as.double(work),  as.integer(ier) )


  ier <- result[[11]]
  if (ier == 1) stop (paste("N = ",n," not valid."))
  if (ier == 2) stop ("A problem with knots detected.")
  if (ier == 3) stop ("Singular coefficient matrix detected.")
  if (ier == 4) stop (paste("NDERIV = ",nderiv," not valid."))
  if (ier == 5) stop (paste("NORDER = ",norder," not valid."))
  if (ier == 6) stop ("X values not strictly increasing.")

  dy  <- matrix(result[[9]],narg,nvar)
  return(dy)
}

  
plot.smooth.Pspline <- function(splobj, ...) {
  if (is.vector(splobj$ysmth) | dim(splobj$ysmth)[[2]] == 1)  
        plot (splobj$x, splobj$ysmth, ...) else
     matplot (splobj$x, splobj$ysmth, ...)
}

lines.smooth.Pspline <- function(splobj, ...) {
  if (is.vector(splobj$ysmth) | dim(splobj$ysmth)[[2]] == 1)  
        lines (splobj$x, splobj$ysmth, ...) else
     matlines (splobj$x, splobj$ysmth, ...)
}

print.smooth.Pspline <- function(x, ...)
{
        if(!is.null(cl <- x$call)) {
                cat("Call:\n")
                dput(cl)
        }
        cat("\nSmoothing Parameter (Spar):", format(x$spar), "\n")
        cat("Equivalent Degrees of Freedom (Df):", format(x$df), "\n")
        cat("GCV Criterion:", format(x$gcv), "\n")
        cat("CV  Criterion:", format(x$cv),  "\n")
        invisible(x)
}


