A. Overview of the Pspline Module

Functions  smooth.Pspline and  predict.smooth.Pspline  are intended to be
similar in functionality and use to standard Splus functions  smooth.spline
and  predict.smooth.spline.  The main difference is that the use can specify
the order m of the derivative that is penalized in the penalized least
squares criterion (using Latex notation)

n^{-1} \sum_i^n [y_i - h(t_i)]^2 + \lambda \int [D^m h(t)]^2 dt

In smooth.spline the order  m  is 2.  

The main motivation for this extension is that one may want to estimate
a derivative of the fitted function  h .  For example, in the 
principal differential analysis technique described in Ramsay (1996) 
a decomposition of curves is developed that estimates a linear differential
operator of order m, and requires derivative estimates through order m.

If derivatives are required, one is generally advised
to let order  m  be two higher than the highest order derivative that is
required, so that the curvature of the derivative is controlled by the 
penalty parameter  \lambda.  For example, if one wanted the second derivative
or acceleration function, one might well let  m = 4.  See the discussion
of Silverman (1985) for many comments on derivative estimation by smoothing
splines.

A derivative of specified order nderiv of h is computed by function 
predict.smooth.Pspline.  The function is called with an argument array that
may not be the same as that used to smooth the original data.  


B. The O(n) Spline Smoothing Algorithm

The algorithm used by smooth.Pspline is an O(n) algorithm, meaning that
the number of operations is proportional to the number of sampled values.
The algorithm was originally defined by Anselone and Laurent (1967)
and later applied to the cubic spline case (m = 2) by Reinsch (1967,1970),
with subsequent refinements by Hutchison and de Hoog (1985).  

The extension to arbitrary m poses a technical problem of computing the 
integral of the product of two basis splines.  While this can be worked out 
analytically, the expressions become horribly complicated very rapidly as 
m increases, and are also inefficient and inaccurate numerically.  
Instead, in this version these integrals are computed using Gaussian 
quadrature, which, since the product of two splines is piecewise a polynomial, 
is exact except for rounding error.  The algorithm also checks for equally 
spaced arguments, in which case the integral needs calculating only once 
instead of n times.

Other applications of the Anselone-Laurent algorithm are described in
a preprint by Heckman and Ramsay (1996).  

There are a few minor differences between smooth.Pspline and smooth.spline
even when  m = 2   that the user should be aware of:

--  smooth.Pspline computes a natural spline while smooth.spline computes
    a B-spline.  The difference between the two will only be apparent at
    the leading and trailing intervals and will generally be small.

--  while smooth.spline permits duplicate argument values,  smooth.Pspline
    does not.  If the data do contain duplicate values for t_i, these
    should be eliminated, the corresponding y_i values averaged,
    and the corresponding weight array w set equal to the number of duplicate
    values.  

--  the method for choosing the smoothing parameter value, called  spar  in
    the argument list, has more options, and is controlled by the  method
    parameter.  Options are:
    method = 1 (default) use the value of spar supplied in the call
    method = 2 adjust  spar  so that the effective degrees of freedom are
       approximately equal to the value of the  df  parameter in the call  
    method = 3 choose  spar  so as to minimize the generalized cross-validation
       criterion (GCV)
    method = 4 choose  spar  so as to minimize the ordinary cross-validation
       criterion (CV)
  
--  minimizing either the GCV or the CV criterion can potentially be a tricky
    business because there are often multiple local minima as well as a global
    minimum.  Functions smooth.spline and smooth.Pspline both use a well-known
    safe-guarded minimization technique referred to as Brent's method
    in Press, et al (1992).  They differ, however, in that smooth.spline
    starts this minimization technique off with an internally generated value,
    whereas smooth.Pspline starts the technique off with the value of the spar
    parameter if this is positive, or with the internal value used in
    smooth.spline otherwise.  This means that the user can conduct a crude
    preliminary grid search to identify the region of the global minimum
    using method=1 before refining this estimate with method=3 or 4.
        
The argument list for smooth.Pspline is very similar to that for smooth.spline,
except that the order of the derivative to be penalized is optionally
specified, and that the type of smoothing parameter estimate is determined
by the method argument as indicated above.  The function heading is as follows:

  smooth.Pspline <- function(x, y, w=rep(1, length(x)), norder=2, 
                           df=norder+2, spar=0, method=1) 

Similarly, the function heading for  predict.smooth.Pspline is

  predict.smooth.Pspline <- function(splobj, xarg, nderiv = 0) 
  

C.  Installation:

Consult the file INSTALL.TXT.

D. Acknowledgements

Preparation of this module was greatly aided by Trevor Hastie, who made available
the code for smooth.spline.  The decision of keep the profile of smooth.Pspline
as close to that of smooth.spline reflects in part my admiration for the
quality of his software.  Nancy Heckman and Bernard Silverman also provided
valuable assistance and advice.

Needless to say, the deficiencies of this software are entirely my own
responsibility.  Problems and comments should be directed to:

Jim Ramsay
Dept. of Psychology
1205 Dr. Penfield Ave.
Montreal, Quebec, Canada
H3A 1B1

email:  ramsay@psych.mcgill.ca
tel:    (514) 398-6123
fax:    (514) 398-4896


References:

Anselone, P. M. and Laurent, P. J. (1967) A general method for the 
    construction of interpolating or smoothing spline-functions.
    Numerische Mathematik, 12, 66-82.

Heckman, N. and Ramsay, J. O. (1996) Some general theory for spline Smoothing,
    McGill University, Unpublished manuscript.

Hutchison, M. F. and de Hoog, F. R. (1985) Smoothing noisy data with 
    spline functions. Numerische Mathematik, 47, 99-106.

Press, W. H., Teukolsky, S. A., Vertling, W. and Flannery, B. P. (1962)
    Numerical Recipes in Fortran: The Art of Scientific Computing,
    Second Edition.  Cambridge: Cambridge University Press.

Reinsch, C. (1967) Smoothing by spline functions. Numerische Mathematik, 10, 
    177-183.

Reinsch, C. (1970) Smoothing by spline functions II. Numerische Mathematik, 16, 
    451-454.

Statistical Sciences, S-PLUS Programmer's Manual, Version 3.2, Seattle:
    StatSci, a division of MathSoft, Inc.

Silverman, B. W. (1985) Some aspects of the spline smoothing approach
    to non-parametric regression curve fitting, Journal of the Royal
    Statistical Society, Series B, 47, 1-52.
