c   ------------------------------------------------------------------------
c  PSPLINE ...  an O(n) spline smoother with penalty on D^m
c     This version can save intermediate results for reruns with new
c     values of smoothing parameter LAMBDA and can compute the
c     GCV, CV, and DF criteria
c
c  This program sets up the necessary two band-structured matrices,
c    and then either proceeds to a smooth calling SPLCAL
c    if the criterion value is fixed,
c    or calls FMM to optimize the smoothing criterion.  In the latter case
c
c  N        ...  number of data points
c  NVAR     ...  number of sets of function values to be smoothed
c  NORDER   ...  order of derivative to be penalized (max. value = 19)
c  X        ...  array of strictly increasing values of length N
c  W        ...  array of positive weights of length N
c  Y        ...  matrix of values to be smoothed of dimension N by NVAR
c  YHAT     ...  matrix of values of    smooths  of dimension N by NVAR
c  LEV      ...  array of N leverage values
c  GCV      ...  output value of the GCV criterion
c  CV       ...  output value of the CV  criterion
c  DF       ...  output value of the  DF criterion
c  LAMBDA   ...  penalty parameter
c  DFMAX    ...  largest tolerated degrees of freedom
c  WORK     ...  working storage array of length at least
c                  (N-NORDER)*(3*NORDER+2)+N
c
c            NB:  if the data are to be re-smoothed with a different
c            value of LAMBDA, the returned contents of WORK must be
c            left untouched between the first and subsequent calls
c
C  METHOD  ... method for computing the smoothing parameter:
c            1  ...  fixed value of LAMBDA
c            2  ...  fixed value of DF
c            3  ...  LAMBDA optimizes GCV criterion
c            4  ...  LAMBDA optimizes  CV criterion
c  IRERUN  ... if nonzero, assume that a previous call has already
c             computed arrays H and GtWG
c  IER     ...  error return:
c            0 ... no error
c            1 ... N < 2*NORDER + 1
c            2 ... NORDER out of permissible range: [1,10]
c            3 ... NVAR < 1
c            4 ... LAMBDA negative
c            5 ... X not strictly increasing
c            6 ... W contains nonpositive elements
c            -j ... failure in the rational Choleski decomposition routine
c                   LDLTBD because of nonpositive value of diagonal element
c                   at index J in the main loop
c
c  See also subroutine SPLIFIT below that evaluates a spline smoothing
c  function or one of its derivatives at user-specified argument values.
c
      subroutine pspline (n, nvar, norder,
     1                x, w, y, yhat, lev, gcv, cv, df,
     2                lambda, dfmax, work, method, irerun, ier)
      implicit double precision(a-h,o-z)
      parameter (NORDLIM = 10)
      dimension x(n), w(n), y(n,nvar), yhat(n,nvar), work(*),
     1          wk1(121), wk2(121)
      double precision lambda, lev(n)
      logical spcwrd
      data eps /1e-7/, tol/1e-3/
c
c  check arguments
c
      if (n .le. 2*norder + 1) then
        ier = 1
        return
      endif
      if (norder .le. 1 .or. norder .gt. NORDLIM) then
        ier = 2
        return
      endif
      if (nvar .lt. 1) then
        ier = 3
        return
      endif
      if (lambda .lt. 0d0) then
        ier = 4
        return
      endif
c
c  Check for x strictly increasing, and also for x being equally spaced.
c  It might save time if this were done prior to calling SPLGCV.
c
      range = x(n) - x(1)
      delta = range/dble(n-1)
      spcwrd = .true.
      critrn = range*eps
      do i=1,n
        if (w(i) .le. 0d0) then
          ier = 6
        endif
        xi = x(i)
        if (spcwrd .and. i .gt. 1 .and.
     1      dabs(xi - xim1 - delta) .gt. critrn) spcwrd = .false.
        if (i .ge. norder .and. xi .le. x(i-norder+1)) then
          ier = 5
        endif
        xim1 = xi
      end do
      if (ier .ne. 0) return
c
      nmnorder = n - norder
c
c  if this is a re-run with a new lambda value, skip computation of H and
c   GtWG
c
      if (irerun .ne. 0) go to 10
c
c  ********************  call HMAT  **************************
c    This step computes band-structured matrix of inner products of
c    B-spline functions of order NORDER
c
      call hmatfn (n, nmnorder, norder, x, work, wk1, spcwrd)
c
c  ********************  call GTWGMAT  **************************
c    This step computes the cross-product
c    of the N by N-NORDER matrix containing the divided difference
c    coefficients, the cross-product having array W as the metric
c
      call gtwgfn (n, norder, x, w, work, wk1, wk2, spcwrd)
c
c  take action depending on METHOD
c
   10 if (method .eq. 1) then
c
c  ********************  call SPLCAL  **************************
c
        call splcal (n, nvar, norder, x, w, y, yhat, lev,
     1               gcv, cv, df, lambda, work, ier)
      else
c
c  ********************  call FMM  **************************
c
	call fmm   (n, nvar, norder, x, w, y, yhat, lev, gcv, cv, df,
     1              lambda, method, work, tol, ier)
        if (ier .ne. 0) return
	if (method .gt. 2 .and. df .gt. dfmax) then
          df = dfmax
          call fmm (n, nvar, norder, x, w, y, yhat, lev, gcv, cv, df,
     1              lambda, 2,      work, tol, ier)
	endif
      endif
c
      return
      end
c   ------------------------------------------------------------------------
c  SPLCAL ...  an O(n) spline smoother with penalty on D^m
c         called by a driver routine that has already set up the two
c         band-structured matrices required in array WORK
c
c  N       ...  number of data points
c  NVAR    ...  number of sets of values to be smoothed
c  NORDER  ...  order of derivative to be penalized (max. value = 19)
c  X       ...  array of strictly increasing values of length N
c  W       ...  array of positive weights of length N
c  Y       ...  N by NVAR matrix of values to be smoothed
c  YHAT    ...  N by NVAR matrix of values of smooth
c  LEV     ...  array of N leverage values
c  GCV     ...  output value of the GCV criterion
c  CV      ...  output value of the CV criterion
c  DF      ...  output value of the DF criterion
c  LAMBDA  ...  penalty parameter
c  WORK    ...  working storage array of length at least
c                  3*N*(NORDER+1)
c
c            NB:  if the data are to be re-smoothed with a different
c            value of LAMBDA, the returned contents of WORK must be
c            left untouched between the first and subsequent calls
c
c  IER     ...  error return:
c            0 ... no error
c            1 ... N < 2
c            2 ... NORDER out of permissible range
c            3 ... X not strictly increasing
c            4 ... LAMBDA negative
c            -j ... failure in the rational Choleski decomposition routine
c                   LDLTBD because of nonpositive value of diagonal element
c                   at index J in the main loop
c
      subroutine splcal (n, nvar, norder, x, w, y, yhat, lev,
     1                   gcv, cv, df, lambda, work, ier)
      implicit double precision(a-h,o-z)
      dimension x(n), w(n), y(n,nvar), yhat(n,nvar), work(*),
     1          wk1(400), wk2(400)
      double precision lambda, lev(n)
c
c  set up offset values for storing information in array WORK
c
      nmnorder = n - norder
      nband    = norder + 1
      nsize    = nmnorder*nband
      igoffset = nmnorder*norder
      iboffset = igoffset + nsize
      iqoffset = iboffset + nsize
      iyoffset = iqoffset + nsize
      nworksiz = nmnorder*(4*norder + 3) + n
c
c    The next result is a band-structured matrix of order N - NORDER
c    and bandwidth  NORDER + 1 of the form  B = H + LAMDA * G^t Wd G
c    where Wd = diag(W) and G contains the divided difference coeffs.
c
      m = 0
      do j=1,norder
        do i=1,nmnorder
          m = m + 1
           work(iboffset+m) = work(m) + lambda*work(igoffset+m)
        end do
      end do
      do i=1,nmnorder
        m = m + 1
        work(iboffset+m) = lambda*work(igoffset+m)
      end do
c
c  ********************  call LBANDMAT  **************************
c    This step computes the rational Choleski decomposition of
c    the above band-structured matrix B
c
      call ldltbdspl (nmnorder, nband, work(iboffset+1), ier)
      if (ier .ne. 0) return
c
c  ++++++++++++++  loop through values to be smoothed  +++++++++++
c
      do ivar=1,nvar
c
c  ********************  call GDIFFFUN  **************************
c    This step computes the divided difference values GtY of Y
c
        do i=1,n
          work(iyoffset+i) = y(i,ivar)
        end do
        call gdifffn (n, norder, x, work(iyoffset+1), wk1, wk2)
c
c  ********************  call SOLVBD  **************************
c    This step solves the equation  BC = GtY.  C replaces GtY
c
        call solvbdspl (nmnorder, nband, 1, work(iboffset+1),
     1                  work(iyoffset+1), ier)
        if (ier .ne. 0) return
c
c  ********************  call GCVEC  **************************
c    This step updates original vector Y to Y - LAMBDA * G X
c
        do i=1,n
          yhat(i,ivar) = y(i,ivar)
        end do
        call gcfn (n, norder, x, w, work(iyoffset+1), yhat(1,ivar),
     1             lambda, wk1, wk2)
      end do
c
c  +++++++++++++++++  end of loop through functions to be smoothed  ++++++
c
c  Compute band-structured portion of inverse of matrix B
c  This replaces B
c
      call bdinvspl (nmnorder, norder, work(iboffset+1), ier)
c     write (*,'(a)') ' B:'
c     do i=1,nmnorder
c       write (*,'(i3,5e12.4)') i,
c    1        (work((j-1)*nmnorder+i+iboffset),j=1,nband)
c     end do
c
c  compute trace of hat matrix, SSE, CV, and GCV criteria
c
      xn    = n
      trace = 0d0
      sse   = 0d0
      cv    = 0d0
      do i=1,n
        sum = 0d0
        ldn = max(0, i-nmnorder)
        lup = min(norder, i-1)
        ml  = nmnorder*ldn
        do l=ldn,lup
          sum = sum + work(ml+i-l+iqoffset)**2*work(i-l+iboffset)
          ml = ml + nmnorder
        end do
        ml1 = nmnorder*ldn
        do l1=ldn,lup-1
          fac = work(ml1+i-l1+iqoffset)
          ml2 = nmnorder*(l1+1)
          do l2=l1+1,lup
            sum = sum +
     1        2d0*fac*work(ml2+i-l2+iqoffset)*
     2        work((l2-l1)*nmnorder+i-l1+iboffset)
            ml2 = ml2 + nmnorder
          end do
          ml1 = ml1 + nmnorder
        end do
        sum = sum*lambda*w(i)
        lev(i) = 1d0 - sum
        trace = trace + sum
        do ivar=1,nvar
          res = (y(i,ivar) - yhat(i,ivar))/w(i)
          sse = sse + res*res
          cv  = cv  + (res/sum)**2
        end do
      end do
c
      gcv = (sse/xn)/(dble(nvar)*trace/xn)**2
      cv  = cv/xn
      df  = xn - trace
c
      return
      end
c  ----------------------------------------------------------------------
c  HMATFN ... computes matrix of inner products of Bspline functions
c
c  X      ...  strictly ascending sequence of values
c  N      ...  length of the sequence
c  NMO    ...  number of rows of band-structured matrix containing inner
c               products:  NMO = N - NORDER.  Number of columns = NORDER
c  NORDER ... order of differential operator in spline penalty term
c             permissible values:  3  ...  4
c  H      ...  band structured matrix with NORDER columns containing nonzero
c             inner products
c  OUTIP  ...  scratch array of length NORDER*(NORDER+1)/2
c  SPCWRD ...  logical flag indicating that X values are equally spaced
c
      subroutine hmatfn (n, nmo, norder, x, h, outip, spcwrd)
      implicit double precision(a-h,o-z)
      dimension x(*), h(nmo,*), outip(*)
      logical spcwrd
c
c  clear array
c
      do i=1,nmo
        do j=1,norder
          h(i,j) = 0d0
        end do
      end do
c
c  ----------------------------  order 1  --------------------------------
c
      if (norder .eq. 1) then
        if (spcwrd) then
          delta = x(2)-x(1)
          do i=1,n-1
            h(i,1) = delta
          end do
        else
          do i=1,n-1
            h(i,1) = x(i+1)-x(i)
          end do
        endif
        return
      endif
c
c  ----------------------------  order 2  --------------------------------
c
      if (norder .eq. 2) then
        if (spcwrd) then
          hi1 = (x(3) - x(1))/3d0
          hi2 = (x(2) - x(1))/6d0
          do i=1,n-2
            h(i,1) = hi1
            if (i .eq. 1) then
              h(i,2) = 0d0
            else
              h(i,2) = hi2
            endif
          end do
        else
          do i=1,n-2
            h(i,1) = (x(i+2)-x(i))/3d0
            if (i .eq. 1) then
              h(i,2) = 0d0
            else
              h(i,2) = (x(i+1) - x(i))/6d0
            endif
          end do
        endif
        return
      endif
c
c  ----------------------------  order 3 or up  ---------------------------
c
      if (norder .gt. 2) then
        nmnorder = n - norder
        if (spcwrd) then
          call splipfn (n, x, norder+1, norder, outip, ier)
          if (ier .ne. 0) then
            ier = ier + 10
            return
          endif
          do i=1,n-1
            m = 0
            do j=1,norder
              imjp1 = i - j + 1
              nmnorderpj = nmnorder + j
              do k=j,norder
                m = m + 1
                kmjp1 = k - j + 1
                if (i .gt. k-1 .and. i .lt. nmnorderpj) then
                  h(imjp1,kmjp1) = h(imjp1,kmjp1) + outip(m)
                endif
              end do
            end do
          end do
        else
          do i=1,n-1
            call splipfn (n, x, i, norder, outip, ier)
            if (ier .ne. 0) then
              ier = ier + 10
              return
            endif
            m = 0
            do j=1,norder
              imjp1 = i - j + 1
              nmnorderpj = nmnorder + j
              do k=j,norder
                m = m + 1
                kmjp1 = k - j + 1
                if (i .gt. k-1 .and. i .lt. nmnorderpj) then
                  h(imjp1,kmjp1) = h(imjp1,kmjp1) + outip(m)
                endif
              end do
            end do
          end do
        endif
      endif
c
      return
      end
c  ----------------------------------------------------------------------
c  GTWGFN ... computes cross-product of divided difference coefficients
c               with respect to x and weights w
c               That is, computes G'WG where G is N by N-NORDER differencing
c               matrix and W is a diagonal weight matrix
c
c  N  ...  length of the sequence
c  NORDER ... order of differential operator in spline penalty term
c             permissible values:  1  ...  4
c  M ...   N - NORDER
c  X  ...  strictly ascending sequence of N values
c  W  ...  N positive weights
c  GTWG  ...  resulting N - NORDER by NORDER + 1 matrix in band structured mode
c  WK    ...  working array of length NORDER
c  C     ...  working array of length NORDER*NORDER
c  SPCWRD ..  logical flag indicating equal spacing of the X values
c
      subroutine gtwgfn (n, norder, x, w, work, wk, c, spcwrd)
      implicit double precision(a-h,o-z)
      dimension x(*), w(*), work(*), c(20,*), wk(*)
      logical spcwrd
c
      nordp1 = norder + 1
      nmnorder = n - norder
      nsize    = nmnorder*nordp1
      igoffset = nmnorder*norder
      iboffset = igoffset + nsize
      iqoffset = iboffset + nsize
c
      if (spcwrd) then
        call divdifffn (nordp1, x(1), c(1,1), wk)
        do i=1,nmnorder
          mj = i
          do j=1,nordp1
            work(iqoffset+mj) = c(j,1)
            mj = mj + nmnorder
          end do
          mj = i
          do j=1,min(i,nordp1)
            sum = 0d0
            do l=1,norder+2-j
              sum = sum + c(l,1)*c(l+j-1,1)*w(i+l-1)
            end do
            work(igoffset+mj) = sum
            mj = mj + nmnorder
          end do
        end do
      else
        do i=1,nmnorder
          call divdifffn (nordp1, x(i), c(1,1), wk)
          mj = i
          do j=1,nordp1
            work(iqoffset+mj) = c(j,1)
            mj = mj + nmnorder
          end do
          mj = i
          do j=1,min(i,nordp1)
            sum = 0d0
            do l=1,norder+2-j
              sum = sum + c(l,1)*c(l+j-1,j)*w(i+l-1)
            end do
            work(igoffset+mj) = sum
            mj = mj + nmnorder
          end do
          do l=1,nordp1
            do j=1,norder
               c(l,norder+2-j) = c(l,nordp1-j)
            end do
          end do
        end do
      endif
c
c  clear upper right triangle
c
      mj = nmnorder
      do j=1,norder
        do i=1,j
          work(mj+i+igoffset) = 0d0
        end do
        mj = mj + nmnorder
      end do
c
c     write (*,'(a)') ' H:'
c     do i=1,nmnorder
c       write (*,'(i3,5f12.4)') i,
c    1        (work((j-1)*nmnorder+i),j=1,norder)
c     end do
c     write (*,'(a)') ' GtG:'
c     do i=1,nmnorder
c       write (*,'(i3,5e12.4)') i,
c    1        (work((j-1)*nmnorder+i+igoffset),j=1,nordp1)
c     end do
c     write (*,'(a)') ' Q:'
c     do i=1,nmnorder
c       write (*,'(i3,5f12.4)') i,
c    1        (work((j-1)*nmnorder+i+iqoffset),j=1,nordp1)
c     end do
c
      return
      end
c  ---------------------------------------------------------------------------
c  SPLIP ...  Computes the inner product of order M bspline functions that
c              are nonzero over an interval.
c              For example, if the lower bound is the
c              knot with index 0, and the order is 4,
c              the active spline functions are:
c              B_{i4}, B_{i-1,4}, B_{i-2,4}, and B_{i-3,4}
c              The inner products are turned in array OUTIP in the order
c              (B_{i4},   B_{i4}),      (B_{i4},   B_{i-1,4}),
c              (B_{i4},   B_{i-2,4}),   (B_{i4},   B_{i-3,4}),
c              (B_{i-1,4},B_{i-1,4}),   (B_{i-1,4},B_{i-2,4}),
c              (B_{i-1,4},B_{i-3,4}),   (B_{i-2,3},B_{i-2,3}),
c              (B_{i-2,4},B_{i-3,4}),   (B_{i-3,4},B_{i-3,4})
c
c  N  ...  length of knot sequence
c  X  ...  strictly increasing sequence of N knot values
c  INDEX ... index of lower bound of interval, must be < N
c  NORDER ...  order of spline < 20
c  OUTIP ...  array of length NORDER*(NORDER+1)/2 for returning inner products
c  IER ... error return
c
      subroutine splipfn (n, x, index, norder, outip, ier)
      implicit double precision (a-h,o-z)
      dimension x(*), outip(*), quadpt(20), quadwt(20), biatx(20)
      double precision knot(40)
c
      ier = 0
      if (index .lt. 1 .or. index .ge. n) then
        ier = 1
        return
      endif
c
c  generate quadrature points and weights for Gauss-Legendre quadrature
c
      call gaulegfn (norder, x(index), x(index+1), quadpt, quadwt)
c
      do m=1,norder*(norder+1)/2
        outip(m) = 0d0
      end do
c
c  first compute local knot sequence from X
c
      knot(norder)   = x(index)
      knot(norder+1) = x(index+1)
      do i=1,norder-1
        if (index-i .ge. 1) then
          knot(norder-i) = x(index-i)
        else
          knot(norder-i) = x(1)
        endif
        if (index+i+1 .le. n) then
          knot(norder+i+1) = x(index+i+1)
        else
          knot(norder+i+1) = x(n)
        endif
      end do
c
c  now compute the spline values at the quadrature points
c
      do i=1,norder
        call bsplvbfn (knot, norder, quadpt(i), norder, biatx)
        m = 0
        wi = quadwt(i)
        do j=1,norder
          do k=j,norder
            m = m + 1
            outip(m) = outip(m)
     1               + wi*biatx(norder-j+1)*biatx(norder-k+1)
          end do
        end do
      end do
c
      return
      end
c ---------------------------------------------------------------------------
c  GAULEG  ...  computes Gauss-Legendre quadrature points and weights
c    for definite integral with unit kernel
c    Adapted from Numerical Recipes in Fortran, p. 145
c
c  Arguments:
c  N     ...  number of points
c  A     ...  lower limit of integration
c  B     ...  upper limit of integration
c  QUADPT ...  quadrature points or abscissas
c  QUADWT ...  quadrature weights
c
      subroutine gaulegfn (n, a, b, quadpt, quadwt)
      implicit double precision (a-h,o-z)
      parameter (EPS=3.0D-14, PI=3.141592654D0)
      dimension quadpt(n), quadwt(n)
c
      dn = n
      m  = (n + 1)/2
      xm = (a + b)/2d0
      x1 = (b - a)/2d0
      do i=1,m
        z = dcos(PI*(dble(i) - 0.25d0)/(dn + 0.5d0))
   10   p1 = 1d0
        p2 = 0d0
        do j=1,n
          dj = j
          p3 = p2
          p2 = p1
          p1 = ((2d0*dj-1d0)*z*p2 - (dj-1d0)*p3)/dj
        end do
        pp = dn*(z*p1-p2)/(z*z-1d0)
        z1 = z
        z  = z1 - p1/pp
        if (dabs(z - z1) .gt. EPS) go to 10
        quadpt(    i) = xm - x1*z
        quadpt(n+1-i) = xm + x1*z
        quadwt(    i) = 2d0*x1/((1d0-z*z)*pp*pp)
        quadwt(n+1-i) = quadwt(i)
      end do
c
      return
      end
c  --------------------------------------------------------------------------
      subroutine bsplvbfn (t, norder, x, left, biatx)
      implicit double precision(a-h, o-z)
      dimension biatx(*), t(*), deltal(20), deltar(20)
c
      j = 1
      biatx(1) = 1d0
      if (j .ge. norder) return
c
   10 jp1 = j + 1
      deltar(j) = t(left+j) - x
      deltal(j) = x - t(left+1-j)
      saved = 0d0
      do i=1,j
         term = biatx(i)/(deltar(i) + deltal(jp1-i))
         biatx(i) = saved + deltar(i)*term
         saved = deltal(jp1-i)*term
      end do
      biatx(jp1) = saved
      j = jp1
      if (j .lt. norder) go to 10
c
      end
c  -------------------------------------------------------------------
c  LDLTBD ...  computes rational Choleski factorization
c            A = LDL' for a banded matrix A
c
c  N  ... order of matrix
c  K  ... number of off-diagonal bands + 1
c  ABAND  ... N by K matrix ... diagonal in 1st column, and
c         column j contains the N - j + 1 nonzero elements
c         of off-diagonal band j starting in row j
c         On exit, the first column contains the values of D, and
c         remaining columns contain the off-diagonal values of L
c  IER  ...  error return:  0 means no error,
c                           1 means N < 1
c                           2 means K < 1
c                           3 means K > N
c                          -J means zero or negative element for D found on
c                             loop J
c
      subroutine ldltbdspl (n, k, aband, ier)
      double precision aband(n,k), vj, sum
c
      do j=1,n
        ist = max (1, j - k + 1)
        do i=ist,j-1
          jmi = j - i
          aband(jmi,k) = aband(j,jmi+1)*aband(i,1)
        end do
        sum = aband(j,1)
        do i=ist,j-1
          jmi = j - i
          sum = sum - aband(j,jmi+1)*aband(jmi,k)
        end do
        vj = sum
        if (vj .le. 0d0) then
          ier = -j
          return
        endif
        aband(j,1) = sum
        iend = min(n,j + k - 1)
        do i=j+1,iend
          sum = aband(i,i-j+1)
          lst = max(1, i - k + 1)
          do l=lst,j-1
            sum = sum - aband(i,i-l+1)*aband(j-l,k)
          end do
          aband(i,i-j+1) = sum/vj
        end do
      end do
c
c  clean up working storage
c
      do i=1,k-1
        aband(i,k) = 0
      end do
c
      return
      end
c  ----------------------------------------------------------------------
c  GDIFFFN ... computes differences for a vector y with respect to x
c               That is, computes G'y where G is N by N-NORDER differencing
c               matrix.
c
c  N  ...  length of the sequence
c  NORDER ... order of differential operator in spline penalty term
c             permissible values:  3  ...  4
c  X  ...  strictly ascending sequence of values
c  Y  ...  sequence to be differenced
c
      subroutine gdifffn (n, norder, x, y, wk, c)
      implicit double precision(a-h,o-z)
      dimension x(*), y(*), wk(*), c(*)
c
      nordp1 = norder + 1
      do i=1,n-norder
        call divdifffn(nordp1, x(i), c, wk)
        sum = 0
          do j=1,nordp1
            sum = sum + y(i+j-1)*c(j)
          end do
        y(i) = sum
      end do
c
      return
      end
c  -------------------------------------------------------------------------
c  DIVDIFFFN  ...  computes divided difference coefficients up to order N
c          for argument value array of length N, N >= 2
c
      subroutine divdifffn (n, x, c, wk)
      implicit double precision (a-h, o-z)
      dimension x(n), c(n), wk(n,n)
c
      if (n .eq. 1) then
        c(1) = 1d0
      endif
c
c  set up coefficients for order 2
c
      nm1 = n - 1
      do i=1,n
        do j=1,nm1
          wk(i,j) = 0d0
        end do
      end do
      do j=1,nm1
        dj1 = x(j+1) - x(j)
        wk(j,j)   = -1d0/dj1
        wk(j+1,j) =  1d0/dj1
      end do
c
c  recurse up to order n - 2
c
      do k=1,n-2
        do j=1,nm1-k
          djk = x(j+k+1) - x(j)
          do i=j,j+k+1
            wk(i,j)  = (wk(i,j+1) - wk(i,j))/djk
          end do
        end do
      end do
c
c  return divided difference coefficients times final difference
c
      do i=1,n
        c(i) = wk(i,1)*djk
      end do
c
      return
      end
c  -------------------------------------------------------------------
c  SOLVBD ...  computes the solution to the equation
c                Ax = y for a symmetric banded matrix A
c            given the rational Choleski factorization  A = LDL'
c
c  N  ... order of matrix
c  K  ... number of off-diagonal bands + 1
c  M  ... number of columns of right side array Y
c  LBAND  ... N by K matrix ... D in 1st column, and
c         column j contains the N - j + 1 nonzero elements of lower triangular
c         Choleski factor L
c         off-diagonal band j starting in row j
c  Y  ... N by M array containing right side.  On return it contains
c         the solutions
c  IER  ...  error return:  0 means no error,
c                           1 means N < 1
c                           2 means K < 1
c                           3 means K > N
c                           J + 10 means zero or negative element for D found on
c                             loop J
c
      subroutine solvbdspl (n, k, m, lband, y, ier)
      double precision lband(n,k), y(n,m), sum
c
c  check arguments
c
      if (n .lt. 1) then
        ier = 1
        return
      endif
      if (k .lt. 1) then
        ier = 2
        return
      endif
      if (k .gt. n) then
        ier = 3
        return
      endif
      if (m .lt. 1) then
        ier = 4
        return
      endif
      do j=1,n
        if (lband(j,1) .le. 0d0) then
          ier = j + 10
          return
        endif
      end do
c
c  Solve  Lu = y
c
      do j=1,n
        ist = max (1, j - k + 1)
        do jrs=1,m
          sum = y(j,jrs)
          do i=ist,j-1
            sum = sum - lband(j,j-i+1)*y(i,jrs)
          end do
          y(j,jrs) = sum
        end do
      end do
c
c  Solve Dv = u
c
      do j=1,n
        do jrs=1,m
          y(j,jrs) = y(j,jrs)/lband(j,1)
        end do
      end do
c
c  Solve  L'x = v
c
      do j=1,n
        jcomp = n - j + 1
        ist = max (1, j - k + 1)
        do jrs=1,m
          sum = y(jcomp,jrs)
          do i=ist,j-1
            icomp = n - i + 1
            sum = sum - lband(icomp,j-i+1)*y(icomp,jrs)
          end do
          y(jcomp,jrs) = sum
        end do
      end do
c
      return
      end
c  ----------------------------------------------------------------------
c  GCFN ... computes  GC where G is N by N-NORDER differencing
c               matrix and C is a vector of length N-NORDER
c
c  N      ...  length of the sequence
c  NORDER ...  order of differential operator in spline penalty term
c              permissible values:  1  ...  4
c  X      ...  strictly ascending sequence of N values
c  W      ...  sequence of N positive weights
c  CVEC   ...  N - NORDER vector
c  y      ...  vector of length N containing values to be smoothed:
c          resulting N vector is Y - lambda*G W C
c
      subroutine gcfn (n, norder, x, w, cvec, y, lambda, wk, c)
      implicit double precision(a-h,o-z)
      dimension x(*), w(*), cvec(*), y(*), c(*), wk(*)
      double precision lambda
c
      nordp1 = norder + 1
      do i=1,n-norder
        factr = lambda*cvec(i)
        call divdifffn (nordp1, x(i), c, wk)
        do l=1,nordp1
          iplm1 = i + l - 1
          y(iplm1) = y(iplm1) - factr*c(l)*w(iplm1)
        end do
      end do
c
      return
      end
!f
c  ---------------------------------------------------------------------
c  BDINV  ...  invert a band-structured matrix of order N and bandwidth
c              M that has been replaced by its rational Choleksi decomp.
c              On completion the inverse overwrites the input matrix X
c
c  N  ...  order of the matrix
c  M  ... number of off-diagonal bands ... M+1 = number of cols of X
c  X  ...  band-structured matrix containing Choleski factors
c  IER  ...  error return
c
      subroutine bdinvspl (n, m, x, ier)
      implicit double precision (a-h,o-z)
      dimension x(n,*)
c
c  check for zero diagonal entries
c
      do i=1,n
         if (x(i,1) .le. 0d0) then
           ier = 10 + i
           return
         endif
      end do
c
      mp1  = m + 1
      ilim = 1
      x(n,1) = 1d0/x(n,1)
      do i=n-1,1,-1
        do l=1,ilim
          sum = 0d0
          ipl = i + l
          do k=1,ilim
            kp1 = k + 1
            ipk = i + k
            if     (k .eq. l) then
              sum = sum - x(ipk,kp1)*x(ipl,1)
            elseif (k .gt. l) then
              sum = sum - x(ipk,kp1)*x(ipk,k-l+1)
            else
              sum = sum - x(ipk,kp1)*x(ipl,l-k+1)
            endif
          end do
          x(l,mp1) = sum
        end do
        sum = 1d0/x(i,1)
        do l=1,ilim
          sum = sum - x(i+l,l+1)*x(l,mp1)
        end do
        x(i,1) = sum
        do l=1,ilim
          x(i+l,l+1) = x(l,mp1)
        end do
        if (ilim .lt. m) ilim = ilim + 1
      end do
c
c  clear upper triangle
c
      do l=1,m
        x(l,mp1) = 0d0
      end do
c
      return
      end
c  ---------------------------------------------------------------------------
      subroutine fmm (n, nvar, norder, xvec, wvec, yvec, yhat, lev,
     1                gcv, cv, df, lambda, method, work, tol, ier)
      implicit double precision(a-h,o-z)
      parameter (XDN = 1d-10, XUP = 3)
      dimension xvec(*), wvec(*), yvec(n,nvar), yhat(n,nvar), work(*)
      double precision lambda, lev(*)
c
      targdf   = df
      nmnorder = n - norder
      igoffset = nmnorder*norder
      t1 = 0d0
      t2 = 0d0
      do i=1,nmnorder
        t1 = t1 + work(i)
        t2 = t2 + work(igoffset+i)
      end do
      ratio = t1/t2
c
      eps = 1d0
   10 eps = eps/2d0
      tol1 = 1d0 + eps
      if (tol1 .gt. 1d0) go to 10
      eps = dsqrt(eps)
c
c  ----------------  initialization of lambda  -------------------------
c
      a = XDN
      b = XUP
c     v = a + 0.382*(b - a)
c     write (*, '(a,e10.3)') ' LAMBDA =', lambda
      if (lambda .le. 0d0) then
        v = 0.75d0
      else
        v = (2d0 + dlog(lambda/ratio)/2.772589d0)/6d0
      endif
      w = v
      x = v
      e = 0d0
      lambda = ratio*dexp(2.772589d0*(6d0*x - 2d0))
c     write (*, '(a,e10.3)') ' LAMBDA =', lambda
c
c  Call 1 to SPLCAL
c
      call splcal (n, nvar, norder, xvec, wvec, yvec, yhat, lev,
     1             gcv, cv, df, lambda, work, ier)
      if (ier .ne. 0) return
c     write (*,'(a, f10.5,e10.3,4f12.4)')
c    1                   ' Call 1:', x, lambda, df, gcv, cv, fu
c
      if (method .eq. 2) fx = (targdf-df)**2
      if (method .eq. 3) fx = gcv
      if (method .eq. 4) fx = cv
      fv = fx
      fw = fx
c
c  --------------------  main loop starts here -------------------------
c
   20 xm   = 0.5d0*(a + b)
      tol1 = eps*dabs(x) + tol/3d0
      tol2 = 2d0*tol1
c
c  check stopping criterion
c
      if (dabs(x - xm) .le. (tol2 - 0.5d0*(b - a))) go to 90
c
c is golden-section necessary?
c
      if (dabs(e) .le. tol1) go to 40
c
c  fit parabola
c
      r = (x - w)*(fx - fv)
      q = (x - v)*(fx - fw)
      p = (x - v)*q - (x - w)*r
      q = 2d0*(q - r)
      if (q .gt. 0d0) p = -p
      q =  dabs(q)
      r = e
      e = d
c
c  is parabola acceptable?
c
      if (dabs(p) .ge. dabs(0.5d0*q*r)) go to 40
      if (p .le. q*(a - x))           go to 40
      if (p .ge. q*(b - x))           go to 40
c
c  a parabolic interpolation step
c
      d = p/q
      u = x + d
c
c  f must not be evaluated too close to a or b
c
      if ((u - a) .lt. tol2 .or. b - u .lt. tol2) then
        if (xm - x .ge. 0d0) then
          d =  tol1
        else
          d = -tol1
        endif
      endif
      go to 50
c
c  a golden-section step
c
   40 if (x .ge. xm) e = a - x
      if (x .lt. xm) e = b - x
      d = 0.382d0*e
c
c  f must not be evaluated too close to x
c
   50 if (dabs(d) .ge. tol1) then
        u = x + d
      else
        if (d .ge. 0d0) then
          u = x + tol1
        else
          u = x - tol1
        endif
      endif
      lambda = ratio*dexp(2.772589d0*(6d0*u - 2d0))
c
c  Call 2 to SPLCAL
c
      call splcal (n, nvar, norder, xvec, wvec, yvec, yhat, lev,
     1             gcv, cv, df, lambda, work, ier)
      if (ier .ne. 0) return
      if (method .eq. 2) fu = (targdf-df)**2
      if (method .eq. 3) fu = gcv
      if (method .eq. 4) fu = cv
c     write (*,'(a, f10.5,e10.3,4f12.4)')
c    1                   ' Call 2:', u, lambda, df, gcv, cv, fu
c
c  update  a, b, v, w, and x
c
      if (fu .gt. fx) go to 60
      if (u  .ge. x)  a = x
      if (u  .lt. x)  b = x
      v  = w
      fv = fw
      w  = x
      fw = fx
      x  = u
      fx = fu
      go to 20
c
   60 if (u  .lt. x)  a = u
      if (u  .ge. x)  b = u
      if (fu .le. fw) go to 70
      if (w  .eq. x)  go to 70
      if (fu .le. fv) go to 80
      if (v  .eq. x)  go to 80
      if (v  .eq. w)  go to 80
      go to 20
c
   70 v  = w
      fv = fw
      w  = u
      fw = fu
      go to 20
c
   80 v  = u
      fv = fu
      go to 20
c
c  -------------------  end of main loop  ------------------------------
c
   90 continue
c
      return
      end

c  -------------------------------------------------------------------------
c  SPLIFIT ...  this subroutine inputs a strictly increasing sequence of
c       arguments X, along with NVAR sets of function values Y, and
c       computes the values of NDERIV derivative values at the NARG argument
c       values in XARG.  The points are interpolated by bsplines of degree  NDEG
c       using the "not-a-knot" condition of de Boor's book, page 55
c
c  Arguments:
c
c  N      ... length of arrays TAU and GTAU
c  NARG   ... number of argument values at which the derivatives are required
c  NVAR   ... number of sets of function values
c  NORDER ... order of B-spline (degree + 1) (minimum value 1)
c  NDERIV ... order of derivative (0 min, NORDER - 2 max)
c  X      ... array of strictly increasing argument values
c  Y      ... matrix of degree N by NVAR of function values at argument values X
c  XARG   ... array length NARG of argument values at which the derivative
c             values are to be computed
c  DY     ... matrix of degree NARG by NVAR of returned derivative values
c  WORK   ... working storage of length (2*NORDER-1)*N + NORDER
c  IER    ...  error return:
c               0 ... no error
c               1 ... inappropriate value of N
c               2 ... problem with knots detected in SPLINT
c               3 ... singularity of coefficient matrix detected in BANFAC
c               4 ... inappropriate value of NDERIV
c               5 ... inappropriate value of NDEG
c               6 ... X not strictly increasing
c
      subroutine splifit (n, narg, nvar, norder, nderiv, x, y, xarg, dy,
     1                    work, ier)
      implicit double precision (a-h,o-z)
      dimension x(n), y(n,*), xarg(narg), dy(narg,*), work(*)
c
      ier    = 0
c
c  check arguments
c
      if (n .le. norder) then
        ier = 1
        return
      endif
      if (nderiv .lt. 0 .or. nderiv .ge. norder) then
        ier = 4
        return
      endif
      if (norder .lt. 1) then
        ier = 5
        return
      endif
      do i=2,n
        if (x(i) .le. x(i-1)) then
          ier = 6
          return
        endif
      end do
c
c  compute offsets for arrays
c
      n2    = 2*n
      iboff = 1
      itoff = n+1
      iqoff = n2+norder+1
c
c  construct knot sequence
c
      do k=1,norder
        work(n+k)  = x(1)
        work(n2+k) = x(n)
      end do
      nhalf = norder/2
      do k=norder+1,n
        work(n+k) = x(k-nhalf)
      end do
c
c  -----------  loop through sets of function values
c
      do ivar=1,nvar
c
c  call SPLINT to get the coefficients for the interpolating bspline
c
        call splint ( x, y(1,ivar), work(itoff), n, norder,
     1                work(iqoff), work(iboff), iflag )
        ier = iflag - 1
        if (ier .ne. 0) return
c
c  go through data computing value of derivative nderiv
c
        do iarg=1,narg
          call dpbvalue(work(itoff), work(iboff), n, norder,
     1                  xarg(iarg), nderiv, dy(iarg,ivar))
        end do
      end do
c
      return
      end
c  --------------------------------------------------------------------------
      subroutine splint ( x, y, t, n, k, q, bcoef, iflag )
c
      double precision bcoef(n), y(n), q(*), t(*), x(n), xi
c
      np1   = n + 1
      km1   = k - 1
      kpkm1 = 2*k - 1
      kpkm2 = 2*km1
      left  = k
      lenq  = n*(k+km1)
      do  i=1,lenq
         q(i) = 0d0
      end do
c
c  ***   loop over i to construct the  n  interpolation equations
c
      do i=1,n
         xi   = x(i)
         ilp1mx = min0(i+k,np1)
c        *** find  left  in the closed interval (i,i+k-1) such that
c                t(left) .le. x(i) .lt. t(left+1)
c        matrix is singular if this is not possible
         left = max0(left,i)
         if (xi .lt. t(left)) then
            iflag = 2
            return
         endif
   10    if (xi .lt. t(left+1))    go to 20
         left = left + 1
         if (left .lt. ilp1mx)       go to 10
         left = left - 1
         if (xi .gt. t(left+1)) then
            iflag = 2
            return
         endif
c
   20    call dpbsplvb ( t, k, 1, xi, left, bcoef )
c
         jj = i-left+1 + (left-k)*(k+km1)
         do j=1,k
            jj = jj+kpkm2
            q(jj) = bcoef(j)
         end do
      end do
c
c     ***obtain factorization of  a  , stored again in  q.
c
      call banfac ( q, kpkm1, n, km1, km1, iflag )
c
      if (iflag .ne. 1) then
         iflag = 3
         return
      endif
c
c     *** solve  a*bcoef = y  by backsubstitution
c
      do i=1,n
         bcoef(i) = y(i)
      end do
c
      call banslv ( q, kpkm1, n, km1, km1, bcoef )
c
      return
      end
c  ---------------------------------------------------------------------
      subroutine banfac ( w, nroww, nrow, nbandl, nbandu, iflag )
c
      integer iflag,nbandl,nbandu,nrow,nroww,   i,ipk,j,jmax,k,kmax
     *                                        ,middle,midmk,nrowm1
      double precision w(nroww,nrow), factor, pivot
c
      iflag = 1
      middle = nbandu + 1
c                         w(middle,.) contains the main diagonal of  a .
      nrowm1 = nrow - 1
c      if (nrowm1)                       999,900,1
      if (nrowm1 .eq. 0)                go to 900
      if (nrowm1 .lt. 0)                go to 999
      if (nbandl .gt. 0)                go to 10
c                a is upper triangular. check that diagonal is nonzero .
      do i=1,nrowm1
         if (w(middle,i) .eq. 0d0)       go to 999
      end do
                                        go to 900
   10 if (nbandu .le. 0) then
c              a is lower triangular. check that diagonal is nonzero and
c                 divide each column by its diagonal .
        do i=1,nrowm1
           pivot = w(middle,i)
           if(pivot .eq. 0d0)              go to 999
           jmax = min0(nbandl, nrow - i)
           do j=1,jmax
              w(middle+j,i) = w(middle+j,i)/pivot
           end do
        end do
        return
      endif
c
c        a  is not just a triangular matrix. construct lu factorization
      do i=1,nrowm1
c                                  w(middle,i)  is pivot for i-th step .
         pivot = w(middle,i)
         if (pivot .eq. 0d0)             go to 999
c                 jmax  is the number of (nonzero) entries in column  i
c                     below the diagonal .
         jmax = min0(nbandl,nrow - i)
c              divide each entry in column  i  below diagonal by pivot .
         do j=1,jmax
            w(middle+j,i) = w(middle+j,i)/pivot
         end do
c                 kmax  is the number of (nonzero) entries in row  i  to
c                     the right of the diagonal .
         kmax = min0(nbandu,nrow - i)
c                  subtract  a(i,i+k)*(i-th column) from (i+k)-th column
c                  (below row  i ) .
         do k=1,kmax
            ipk = i + k
            midmk = middle - k
            factor = w(midmk,ipk)
            do j=1,jmax
               w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)*factor
            end do
         end do
      end do
c                                       check the last diagonal entry .
  900 if (w(middle,nrow) .ne. 0d0)       return
  999 iflag = 2
                                        return
      end
c  ---------------------------------------------------------------------
      subroutine banslv ( w, nroww, nrow, nbandl, nbandu, b )
c
      integer nbandl,nbandu,nrow,nroww,   i,j,jmax,middle,nrowm1
      double precision w(nroww,nrow),b(nrow)
c
      middle = nbandu + 1
      if (nrow .eq. 1)                  go to 20
      nrowm1 = nrow - 1
      if (nbandl .ne. 0) then
c                                 forward pass
c            for i=1,2,...,nrow-1, subtract  right side(i)*(i-th column
c            of  l )  from right side  (below i-th row) .
        do i=1,nrowm1
           jmax = min0(nbandl, nrow-i)
           do j=1,jmax
              b(i+j) = b(i+j) - b(i)*w(middle+j,i)
           end do
        end do
      endif
c                                 backward pass
c            for i=nrow,nrow-1,...,1, divide right side(i) by i-th diag-
c            onal entry of  u, then subtract  right side(i)*(i-th column
c            of  u)  from right side  (above i-th row).
      if (nbandu .le. 0) then
c                                a  is lower triangular .
        do i=1,nrow
           b(i) = b(i)/w(1,i)
        end do
        return
      endif
      i = nrow
   10    b(i) = b(i)/w(middle,i)
         jmax = min0(nbandu,i-1)
         do j=1,jmax
            b(i-j) = b(i-j) - b(i)*w(middle-j,i)
         end do
         i = i - 1
         if (i .gt. 1)                  go to 10
   20 b(1) = b(1)/w(middle,1)
c
      return
      end
c  -------------------------------------------------------------------------
      subroutine dpbsplvb ( t, jhigh, index, x, left, biatx)
c
      parameter (jmax = 20)
      integer index,jhigh,left,   i,j,jp1
      double precision biatx(jhigh), t(*), x, deltal(jmax),
     1     deltar(jmax), saved, term
c
      data j/1/
      save j,deltal,deltar
c
                                        go to (10,20), index
   10 j = 1
      biatx(1) = 1d0
      if (j .ge. jhigh)                 go to 99
c
   20    jp1 = j + 1
         deltar(j) = t(left+j) - x
         deltal(j) = x - t(left+1-j)
         saved = 0d0
         do i=1,j
            term = biatx(i)/(deltar(i) + deltal(jp1-i))
            biatx(i) = saved + deltar(i)*term
            saved = deltal(jp1-i)*term
         end do
         biatx(jp1) = saved
         j = jp1
         if (j .lt. jhigh)              go to 20
c
   99                                   return
      end
c  --------------------------------------------------------------------------
      subroutine dpbvalue ( t, bcoef, n, k, x, jderiv, fofx )
c
      double precision bcoef(n), t(*), x,  aj(20), dl(20), dr(20),
     1  fkmj, fofx
c
      fofx = 0.d0
      if (jderiv .ge. k) return
c
      call dpinterv ( t, n+k, x, i, mflag )
c     if (mflag .ne. 0) return
c  *** if k = 1 (and jderiv = 0), fofx = bcoef(i).
      km1 = k - 1
      if (km1 .eq. 0) then
        fofx = bcoef(i)
        return
      endif
c
c  *** store the k b-spline coefficients relevant for the knot dpinterval
c     (t(i),t(i+1)) in aj(1),...,aj(k) and compute dl(j) = x - t(i+1-j),
c     dr(j) = t(i+j) - x, j=1,...,k-1 . set any of the aj not obtainable
c     from input to zero. set any t.s not obtainable equal to t(1) or
c     to t(n+k) appropriately.
      jcmin = 1
      imk   = i - k
      ip1   = i + 1
      if (imk .lt. 0) then
        jcmin = 1 - imk
        do j=1,i
           dl(j) = x - t(ip1-j)
        end do
        do j=i,km1
           aj(k-j) = 0d0
           dl(j)   = dl(i)
        end do
      else
        do j=1,km1
           dl(j) = x - t(ip1-j)
        end do
      endif
c
      jcmax = k
      nmi   = n - i
      if (nmi .lt. 0) then
        jcmax = k + nmi
        do j=1,jcmax
           dr(j) = t(i+j) - x
        end do
        do j=jcmax,km1
           aj(j+1) = 0d0
           dr(j)   = dr(jcmax)
        end do
      else
        do j=1,km1
           dr(j) = t(i+j) - x
        end do
      endif
c
      do jc=jcmin,jcmax
         aj(jc) = bcoef(imk + jc)
      end do
c
c               *** difference the coefficients  jderiv  times.
      if (jderiv .gt. 0) then
        do j=1,jderiv
           kmj  = k-j
           fkmj = dble(kmj)
           ilo  = kmj
           do jj=1,kmj
             aj(jj) = ((aj(jj+1) - aj(jj))/(dl(ilo) + dr(jj)))*fkmj
             ilo = ilo - 1
           end do
        end do
      endif
c
c  *** compute value at  x  in (t(i),t(i+1)) of jderiv-th derivative,
c     given its relevant b-spline coeffs in aj(1),...,aj(k-jderiv).
      if (jderiv .lt. km1) then
        jdrvp1 = jderiv + 1
        do j=jdrvp1,km1
           kmj = k-j
           ilo = kmj
           do jj=1,kmj
              aj(jj) = (aj(jj+1)*dl(ilo) + aj(jj)*dr(jj))
     1                 /(dl(ilo)+dr(jj))
              ilo = ilo - 1
           end do
        end do
      endif
c
      fofx = aj(1)
c
      return
      end
c  --------------------------------------------------------------------------
      subroutine dpinterv ( xt, lxt, x, left, mflag )
c
      integer left,lxt,mflag,   ihi,ilo,istep,middle
      double precision x,xt(lxt)
c
      data ilo /1/
      save ilo
c
      ihi = ilo + 1
      if (ihi .lt. lxt)                 go to 20
         if (x .ge. xt(lxt))            go to 110
         if (lxt .le. 1)                go to 90
         ilo = lxt - 1
         ihi = lxt
c
   20 if (x .ge. xt(ihi))               go to 40
      if (x .ge. xt(ilo))               go to 100
c
c              **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
      istep = 1
   31    ihi = ilo
         ilo = ihi - istep
         if (ilo .le. 1)                go to 35
         if (x .ge. xt(ilo))            go to 50
         istep = istep*2
                                        go to 31
   35 ilo = 1
      if (x .lt. xt(1))                 go to 90
                                        go to 50
c              **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
   40 istep = 1
   41    ilo = ihi
         ihi = ilo + istep
         if (ihi .ge. lxt)              go to 45
         if (x .lt. xt(ihi))            go to 50
         istep = istep*2
                                        go to 41
   45 if (x .ge. xt(lxt))               go to 110
      ihi = lxt
c
c           **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the dpinterval.
   50 middle = (ilo + ihi)/2
      if (middle .eq. ilo)              go to 100
c     note. it is assumed that middle = ilo in case ihi = ilo+1 .
      if (x .lt. xt(middle))            go to 53
         ilo = middle
                                        go to 50
   53    ihi = middle
                                        go to 50
c**** set output and return.
   90 mflag = -1
      left = 1
                                        return
  100 mflag = 0
      left = ilo
                                        return
  110 mflag = 1
	  if (x .eq. xt(lxt)) mflag = 0
      left = lxt
  111 if (left .eq. 1)                  return
	  left = left - 1
	  if (xt(left) .lt. xt(lxt))        return
                                            go to 111
      end
