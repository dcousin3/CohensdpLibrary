function klprimesmixtureidf( p, n, obsr, obsd, delta, maxitr, ier )

!-----------------------------------------------------------------------
!
!     Calculates the quantile oft a random variable distributed
!     according to the L' distribution with parameters 
!           n        degree of freedom, and
!           obsd     centrality parameter itself a function of rho 
!     where rho is distributed according to the K' distribution with parameters
!           n-2, n-1 degrees of freedom and
!           obsr     centrality parameter
!
!     X     - Input . delta value where to compute the pdf     - Real
!     N     - Input . Sample size                     (N >  0) - Integer
!     OBSR  - Input . for eccentricity parameter of K'         - Real
!     OBSD  - Input . for eccentricity parameter of L'         - Real
!     DELTA - Input . Maximum absolute error required on       - Real
!                     kprimecdf (stopping criteria)
!                     (eps < DELTA < 1 where eps is machine
!                     epsilon; see parameter statement below)
!     MAXITR- Input . Maximum number of iterations            - Integer
!     IER   - Output. Unreturned                              - Integer
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

    !  Function
    !  --------
    real(kind=8) :: klprimesmixtureidf

    !  Arguments
    !  ---------
    real(kind=8), intent(in)  :: p, n, obsr, obsd, delta
    integer,      intent(in)  :: maxitr
    integer,      intent(out) :: ier

    !  Local declarations
    !  ------------------
    real(kind=8)              :: k, ncL, dlL, scaled
    real(kind=8), external    :: kprimelprimeconvolutioncdf, mydtrinv

    dlL = 2 * (n- 1._8) / (1+obsr**2)
    scaled = sqrt(n / (2 * (1._8- obsr)) )
    ncL = scaled * obsd 
    k = sqrt(2._8 / dlL) * Exp(dlgama((dlL+1)/2._8) - dlgama(dlL/2._8) )

    klprimesmixtureidf = mydtrinv( func, p,       &
                .FALSE., .FALSE., -1.0e6_8, 1.0e6_8,      &  
                ncL * k / scaled,                         & ! mean of unperturbated L' distribution
                sqrt( ncL**2 * (1-k**2) + 1 )/scaled ,    & ! std  of unperturbated L' distribution
                +1.0e-4_8, 20,                            &
                ier                                       &
            ) 

contains
    function func(x, iok)
        real(kind=8), intent(in) :: x
        integer,      intent(out):: iok
        real(kind=8), external   :: kprimelprimeconvolutioncdf
        real(kind=8) :: func
        func = kprimelprimeconvolutioncdf(x, n, obsr, obsd, delta, maxitr, iok)
    end function func

end function klprimesmixtureidf


