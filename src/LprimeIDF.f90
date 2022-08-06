FUNCTION lprimeidf( p, q, a1, TOL, MAXITER, ier )
    !-----------------------------------------------------------------------
    !     Calculates the quantile of a random variable distributed
    !     according to the Lambda' distribution with Q degrees of freedom,
    !     A centrality parameter, that is the X such that F( X | q, a) = p
    !
    !     P     - Input . Value of the quantile          (0<P<1)  - Real
    !     Q     - Input . Degrees of freedom             (Q >  0) - Real
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
    !       LPRIMECDF
    !     Fortran functions called:
    !       ABS  LOG  SQRT  DLGAMA
    !
    !*********************************************************************************************!
    !**                                                                                         **!
    !** This function was added by Denis Cousineau, 28 november 2020.                           **!
    !** It is just a wrapper to the generic function dtrinv.f90 made by J. Poitevineau.         **!
    !**                                                                                         **!
    !*********************************************************************************************!

    IMPLICIT NONE
    INTEGER, PARAMETER        :: PR=KIND(1.0D0)

    !  Function
    REAL(PR) :: lprimeidf

    !  Arguments
    REAL(PR), INTENT(in)  :: p, q, a1, TOL
    INTEGER,  INTENT(in)  :: MAXITER
    INTEGER,  INTENT(out) :: ier

    !  Local declarations
    REAL(PR) :: k
    REAL(PR), EXTERNAL :: lprimecdf, dtrinv

    ier = 0

    k = Exp(dlgama((q+1)/2.) - dlgama(q/2.) ) * sqrt(2/q)

    lprimeidf = dtrinv( func, p,                       &
                .FALSE., .FALSE., -1.0D6, 1.0D6,       &
                a1 * k,                                &
                sqrt(a1**2 *(1-k**2) + 1 ),            &
                +1.0D-6, 1000,                         &
                ier                                    &
            )

CONTAINS
    FUNCTION func(x, iok)
        REAL(PR), INTENT(in) :: x
        INTEGER,  INTENT(out):: iok
        REAL(PR), EXTERNAL   :: lprimecdf
        REAL(PR) :: func
        iok = 0
        func = lprimecdf(x, q, a1, TOL, MAXITER, iok)
    END FUNCTION func

END FUNCTION lprimeidf


