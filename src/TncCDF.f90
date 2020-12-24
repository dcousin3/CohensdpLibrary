function tnccdf( x, q, a, delta, maxitr, ier )

!-----------------------------------------------------------------------
!
!     Computes the probability that a random variable distributed
!     according to the noncentral t distribution with Q degrees of
!     freedom, A centrality parameter, is less than or equal to X
!     Note: P( t'q (x) < a ) = P( L'q (a) > x )
!
!     X     - Input . Value of the variable                   - Real
!     Q     - Input . Degrees of freedom             (Q >  0) - Real
!     A     - Input . Eccentricity parameter                  - Real
!     DELTA - Input . Maximum absolute error required on      - Real
!                     tnccdf (stopping criterion)
!                     (eps < delta < 1 where eps is machine
!                     epsilon; see parameter statement below)
!     MAXITR- Input . Maximum number of iterations            - Integer
!     IER   - Output. Return code :                           - Integer
!                     0 = normal
!                     1 = invalid input argument
!                         (then tnccdf = one)
!                     2 = maximum number of iterations reached
!                         (then tnccdf = value at last iteration,
!                                    or one)
!                     3 = required accuracy cannot be reached
!                         (then tnccdf = value at last iteration)
!                     4 = error in auxiliary function (chi2cdf or tcdf)
!                     7 = result out of limits (i.e. <0 or >1)
!
!     External functions called:
!       LPRIMECDF
!
!-----------------------------------------------------------------------

   implicit none

   !  Function
   !  --------

   real(kind=8) :: tnccdf

   !  Arguments
   !  ---------

   real(kind=8), intent(in) :: x, q, a, delta
   integer, intent(in) :: maxitr
   integer, intent(out) :: ier

   !  local declarations
   !  ------------------

   real(kind=8), external :: lprimecdf

   real(kind=8), parameter :: one=1.0_8

   !-----------------------------------------------------------------

   tnccdf = one - lprimecdf( a, q, x, delta, maxitr, ier )

end function tnccdf
