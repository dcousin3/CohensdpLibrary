function kprimeidf( p, q, r, a1, TOL, MAXITER, ier )

!-----------------------------------------------------------------------
!
!     Calculates the probability that a random variable distributed
!     according to the K' distribution with Q and R degrees of
!     freedom, A1 centrality parameter, is less than or equal to X
!
!     P     - Input . p value of the desired quantile (0<P<1) - Real
!     Q     - Input . First degrees of freedom        (Q > 0) - Real
!     R     - Input . Second   "    "     "           (R > 0) - Real
!     A1    - Input . Eccentricity parameter                  - Real
!     TOL   - Input . Maximum absolute error required on      - Real
!                     kprimecdf (stopping criteria)
!                     (eps < TOL < 1 where eps is machine
!                     epsilon; see parameter statement below)
!     MAXITER- Input . Maximum number of iterations            - Integer
!     IER   - Output. Return code :                           - Integer
!                     0 = normal
!                    -1 = no more evolution of the sum but
!                         required accuracy not reached yet
!                         (then kprimecdf = value at last iteration)
!                     1 = invalid input argument
!                         (then kprimecdf = zero)
!                     2 = maximum number of iterations reached
!                         (then kprimecdf = value at last iteration)
!                     3 = cannot be computed
!                         (then kprimecdf = zero)
!                     4 = error in auxiliary function
!                         (betacdf, lprimecdf or tcdf)
!                     5 = result out of limits (i.e. <0 or >1)
!                     7 = 2 + 5 both codes apply
!
!     External functions called:
!       KPRIMECDF
!     Fortran functions called:
!       ABS  LOG  SQRT  DLGAMA
!
!*********************************************************************************************!
!**                                                                                         **!
!** This function was added by Denis Cousineau, 28 november 2020.                           **!
!** It is just a wrapper to the generic function dtrinv.f90 made by J. Poitevineau.         **!
!**                                                                                         **!
!*********************************************************************************************!

    implicit none
    INTEGER, PARAMETER        :: PR=KIND(1.0D0)

    !  Function
    !  --------
    real(PR) :: kprimeidf

    !  Arguments
    !  ---------
    real(PR), intent(in)      :: p, q, r, a1, TOL
    integer,  intent(in)      :: MAXITER
    integer,  intent(out)     :: ier

    !  Local declarations
    !  ------------------
    real(PR)           :: k
    real(PR), external :: kprimecdf, dtrinv

    k = Exp(dlgama((q+1)/2.) + dlgama((r-1)/2.) - dlgama(q/2.) - dlgama(r/2.)) * sqrt(r/q)

    kprimeidf = dtrinv( func, p,                            &
                .FALSE., .FALSE., -1.0D6, +1.0D6,           &
                a1 * k,                                     &
                sqrt( (r/(r-2) - k**2) * a1**2 + r/(r-2) ), &
                TOL, 100,                               &
                ier                                         &
            )

contains
    function func(x, iok)
        real(PR), intent(in) :: x
        integer,      intent(out):: iok
        real(PR), external   :: kprimecdf
        real(PR) :: func
        func = kprimecdf(x, q, r, a1, TOL / 10, MAXITER, iok)
    end function func

end function kprimeidf


