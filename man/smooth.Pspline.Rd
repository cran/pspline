\name{smooth.Pspline}
\alias{smooth.Pspline}
\title{
Fit a Polynomial Smoothing Spline of Arbitrary Order
}
\description{
Returns an object of class \code{smooth.Pspline} which is a natural 
polynomial smooth of the input data of order fixed by the user.
}
\usage{
smooth.Pspline(x, y, w=rep(1, length(x)), norder=2, df=norder + 2, 
               spar=0, method=1)
}
\arguments{
\item{x}{
values of the predictor variable.  These must be strictly increasing,
and there must be at least \code{2*norder + 1} of them.
}
\item{y}{
one or more sets of response variable values.  If there is one
response variable, \code{y} is an array of the same length as \code{x}; if more
than one, then \code{y} is a matrix with \code{length(x)} rows and number of
columns equal to the number of variables.
}
\item{w}{
vector of positive weights for smoothing of the same length as \code{x}.
If measurements at different values of \code{x} have different variances, \code{w}
should be inversely proportional to the variances.  The default is
that all weights are one.
}
\item{norder}{
the order of the spline.  \code{norder = 2} gives the cubic smoothing
spline, and more generally the smoothing function is a piecewise
polynomial of degree \code{2*norder - 1}.  If derivatives are to be
computed from the smoothing using \code{predict.smooth.Pspline}, the order
should be one or two more than the highest order of derivative.
}
\item{df}{
a number which specifies the degrees of freedom = trace(S).  Here S is
the implicit smoothing matrix.  \code{df} controls the amount of smoothing
if \code{method = 2}.
}
\item{spar}{
the usual smoothing parameter for smoothing splines, which is the
coefficient of the integrated squared derivative of order \code{norder}.
\code{spar} controls the amount of smoothing if \code{method = 1}.
}
\item{method}{
the method for controlling the amount of smoothing.  \code{method = 1} uses
the value supplied for \code{spar}.  \code{method = 2} adjusts \code{spar} so that the
degrees of freedom is equal to \code{df}.  \code{method = 3} adjusts \code{spar} so
that the generalized cross-validation criterion is minimized.  `method
= 4' adjusts \code{spar} so that the ordinary cross-validation criterion is
minimized.  If \code{method = 3} or \code{method = 4}, \code{spar} defines the initial
value for the minimization algorithm if positive; otherwise an
internally generated value is used.
}
\value{
an object of class \code{"smooth.Pspline"} is returned, consisting of the fitted
smoothing spline evaluated at the supplied data, some fitting criteria
and constants.  This object contains the information necessary to evaluate
the smoothing spline or one of its derivatives at arbitrary argument
values using \code{predict.smooth.Pspline}.  The components of the returned
list are
}
\item{norder}{
the order of the spline  
}
\item{x}{
values of the predictor variable
}
\item{ysmth}{
a matrix with \code{length(x)} rows, each column of which contains
the smoothed response variable values for the corresponding column of \code{y}.  
}
\item{lev}{
leverage values, which are the diagonal elements of the smoother matrix S.
}
\item{gcv}{
generalized cross-validation criterion value
}
\item{cv}{
ordinary cross-validation criterion value
}
\item{df}{
a number which supplies the degrees of freedom = trace(S) rather than
a smoothing parameter.
}
\item{spar}{
the final smoothing parameter for smoothing splines.  This
is unchanged if \code{method = 1}, but adjusted otherwise.
}
\item{call}{
the call that produced the fit.
}}
\details{
The method produces results similar to function \code{smooth.spline}, but
the smoothing function is a natural smoothing spline rather than a B-spline
smooth, and as a consequence will differ slightly for \code{norder = 2} over the
initial and final intervals.  


The main extension is the possibility of setting the order of
derivative to be penalized, so that derivatives of any order can be
computed using the companion function \code{predict.smooth.Pspline}.  The
algorithm is of order N, meaning that the number of floating point
operations is proportional to the number of values being smoothed.
Note that the argument values must be strictly increasing, a condition
that is not required by \code{smooth.spline}.


Note that the appropriate or minimized value of the smoothing parameter
\code{spar} will depend heavily on the order; the larger the order, the smaller
this parameter will tend to be.
}
\references{
Heckman, N. and Ramsay, J. O. (1996) Spline smoothing with model based
penalties.  McGill University, unpublished manuscript.
}
\seealso{
\code{\link{predict.smooth.Pspline}}, \code{\link{plot.smooth.Pspline}}
}
\keyword{smooth}
% Converted by Sd2Rd version 0.3-1.
