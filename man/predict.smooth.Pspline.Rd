\name{predict.smooth.Pspline}
\alias{predict.smooth.Pspline}
\title{
Smoothing Spline of Arbitrary Order at New Data
}
\description{
Uses an object of class \code{"smooth.Pspline"} to evaluate a polynomial smoothing
spline of arbitrary order or one of its derivatives at new argument values.
}
\usage{
predict.smooth.Pspline(splobj, xarg, nderiv=0)
}
\arguments{
\item{splobj}{
a fitted \code{smooth.Pspline} object.
}
\item{xarg}{
the argument values at which the spline or its derivative is to be
evaluated.
}
\item{deriv}{
the order of the derivative required -- default is 0, the function itself.
}}
\value{
a list with components \code{xarg} and \code{dy}; the \code{xarg} component is identical
to the input \code{xarg} sequence, the \code{dy} component is the evaluated derivative
of order \code{deriv}. 
}
\details{
The method produces results similar to function \code{predict.smooth.spline}, but
the smoothing function is a natural smoothing spline rather than a B-spline
smooth, and the order of the spline can be chosen freely, where order
in this case means the order of the derivative that is penalized.  
\code{smooth.spline} penalizes the second derivative, and consequently only
derivatives or order 0 or 1 are useful, but because \code{smooth.Pspline}
penalizes a derivative of order m, derivatives up to order m-1 are useful.
The general recommendation is to penalize the derivative two beyond the
highest order derivative to be evaluated.


}
\references{
Heckman, N. and Ramsay, J. O. (1996) Spline smoothing with model based
penalties.  McGill University, unpublished manuscript.
}
\seealso{
\code{\link{smooth.Pspline}}
}
\examples{
example(smooth.Pspline)
## smoother line is given by
xx <- seq(4, 25, length=100)
lines(xx, predict(sm.spline(speed, dist, df=5), xx), col = "red")
## add plots of derivatives
lines(xx, 10*predict(sm.spline(speed, dist), xx, 1), col = "blue")
lines(xx, 100*predict(sm.spline(speed, dist), xx, 2), col = "green")
}
\keyword{smooth}
