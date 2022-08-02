function tprimecdf( x, q, a, TOL, MAXITER, ier )

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
!     TOL   - Input . Maximum absolute error required on      - Real
!                     tnccdf (stopping criterion)
!                     (eps < TOL < 1 where eps is machine
!                     epsilon; see parameter statement below)
!     MAXITER- Input . Maximum number of iterations            - Integer
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
   INTEGER, PARAMETER        :: PR=KIND(1.0D0)

   !  Function
   !  --------

   real(PR) :: tprimecdf

   !  Arguments
   !  ---------

   real(PR), intent(in) :: x, q, a, TOL
   integer, intent(in)  :: MAXITER
   integer, intent(out) :: ier

   !  local declarations
   !  ------------------

   real(PR), external :: lprimecdf

   real(PR), parameter :: one=1.0D0

   !-----------------------------------------------------------------

   tprimecdf = one - lprimecdf( a, q, x, TOL, MAXITER, ier )

end function tprimecdf
