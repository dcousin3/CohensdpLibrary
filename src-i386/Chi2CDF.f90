function chi2cdf( x, p, plimit, ier )

!-----------------------------------------------------------------------
!
!     Computes the probability that a random variable distributed
!     according to the chi-square distribution with P degrees of
!     freedom, is less than or equal to X.
!
!     X     - Input . Value of the variable         (X >= 0)  - Real
!     P     - Input . Degrees of freedom            (P >  0)  - Real
!     PLIMIT- Input . Use Wilson and Hilferty's     (P >  0)  - Real
!                     approximation when P > PLIMIT
!                     (should be at least 1000)
!     IER   - Output. Return code :                           - Integer
!                     0 = normal
!                     1 = invalid input argument
!                         (then CHI2CDF = -1.)
!                     2 = maximum number of iterations reached
!                         (then CHI2CDF = value at last iteration)
!
!     External functions called:
!       NCDF standard normal distribution
!       DLGAMA (Log(Gamma(.)) if not in FORTRAN library
!
!     Fortran functions called:
!       ABS  EXP  LOG  SQRT  (and DLGAMA if available)
!
!-----------------------------------------------------------------------

   implicit none
   INTEGER, PARAMETER        :: PR=KIND(1.0D0)

   !  Function
   !  --------

   real(PR) :: chi2cdf

   !  Arguments
   !  ---------

   real(PR), intent(in) :: x, p, plimit
   integer, intent(out) :: ier

   !  Local declarations
   !  ------------------

!! real(PR), external :: dlgama
   real(PR), external :: ncdf

   real(PR), parameter :: zero=0.0D0, half=0.5D0, one=1.0D0
   real(PR), parameter :: two=one+one, three=two+one, nine=three*three
   real(PR), parameter :: xlimit=1.0D-40
   real(PR), parameter :: relerr=1.0D-14
   real(PR), parameter :: tinyr=1.0D-307, explower=-706.893D0
   real(PR), parameter :: fpmin=1.0D-300
   integer, parameter :: itmax=2000

   !     relerr = relative error
   !        These constants are machine dependent:
   !     fpmin = a real nearby the smallest positive real
   !     tinyr = the smallest positive real
   !     explower = minimum valid argument for the exponential function
   !             (i.e. log(tinyr))

   real(PR) :: a, an, ap, b, c, d, del, h, y
   integer :: i, n
   logical :: flag

   !--------------------------------------------------------------------

   !  Test for valid input arguments

   if ( x < zero .or. p <= zero .or. plimit <= zero ) then
      ier = 1
      chi2cdf = -one
      return
   end if

   ier = 0
   if ( x < xlimit ) then
      chi2cdf = zero
      return
   end if

   if ( p <= plimit ) then
      y = x * half
      a = p * half

      if ( y < a+one ) then
         flag = .false.
         ap   = a
         h    = one/a
         del  = h
         do n = 1, itmax
            ap  = ap + one
            del = del*y/ap
            h   = h + del
            if ( abs(del) < abs(h)*relerr ) goto 10
         end do
      else
         flag = .true.
         b    = y + one - a
         c    = one/fpmin
         d    = one/b
         h    = d
         do i = 1, itmax
            an = -i*(i-a)
            b  = b + two
            d  = an*d + b
            if ( abs(d) < fpmin ) d = fpmin
            c = b + an/c
            if ( abs(c) < fpmin ) c = fpmin
            d   = one/d
            del = d*c
            h   = h*del
            if ( abs(del-one) < relerr ) goto 10
         end do
      end if
      !  Maximum number of iterations reached
      ier = 2

      10 continue
      a = - y + a*log(y) - dlgama(a)
      if ( h >= tinyr ) then
         a = log(h) + a
         if ( a >= explower ) then
            chi2cdf = exp(a)
         else
            chi2cdf = zero
         end if
      else
         chi2cdf = zero
      end if
      if ( flag ) chi2cdf = one - chi2cdf

   else
      !  Use Wilson & Hilferty's approximation when great df's
      !  Compute proba(U<=H) where U follows a normal(0,1) distr.

      b = two/(nine*p)
      a = log(x/p)/three
      if ( a >= explower ) then
         h = (exp(a)-one+b) / sqrt(b)
      else
         h = (b-one) / sqrt(b)
      end if
      chi2cdf = ncdf(h)
   end if

end function chi2cdf
