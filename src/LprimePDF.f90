function lprimepdf( x, q, a1, delta, maxitr, ier )
    
!-----------------------------------------------------------------------
!
!     Calculates the density that a random variable distributed
!     according to the L' distribution with Q degrees of
!     freedom, A1 centrality parameter, is equal to X
!
!     P     - Input . p value of the desired quantile (0<P<1) - Real
!     Q     - Input . First degrees of freedom       (Q >  0) - Real
!     A1    - Input . Eccentricity parameter                  - Real
!     DELTA - Input . Maximum absolute error required on      - Real
!                     kprimecdf (stopping criteria)
!                     (eps < DELTA < 1 where eps is machine
!                     epsilon; see parameter statement below)
!     MAXITR- Input . Maximum number of iterations            - Integer
!     IER   - Output. unreturned...                           - Integer
!
!     External functions called:
!       LPRIMECDF
!     Fortran functions called:
!       ABS    MAX
!
!*********************************************************************************************!
!**                                                                                         **!
!** This function was added by Denis Cousineau, 28 november 2020.                           **!
!** It is just a wrapper to the generic function dfridr from Numerical Reciepes.            **!
!** Press, Teukolsky, Vetterling, Flannery (1992)                                           **!
!**                                                                                         **!
!*********************************************************************************************!

    implicit none

    !  Function
    !  --------
    real(kind=8) :: lprimepdf

    !  Arguments
    !  ---------
    real(kind=8), intent(in)  :: x, q, a1, delta
    integer,      intent(in)  :: maxitr
    real(kind=8), intent(out) :: ier

    !  Local declarations
    !  ------------------
    real(kind=8), external  :: lprimecdf
    real(kind=8)            :: rer  ! real-valued error 

    ier = 0
    lprimepdf = dfridr( func, x, 0.1_8, rer )

contains
    function func( x )
        real(kind=8), intent(in) :: x
        real(kind=8), external   :: lprimecdf
        real(kind=8)  :: func
        integer       :: iok
        func = lprimecdf(x, q, a1, delta, maxitr, iok)
    end function func

    function dfridr(func, x, h, rer )
        ! Reference: Press, Teukolsky, Vetterling, Flannery (1992) Numerical Receipes in fortran 77 (vol. 1)
        real(kind=8)              :: dfridr
        real(kind=8), external    :: func
        real(kind=8), intent(in)  :: x, h
        real(kind=8), intent(out) :: rer
        real(kind=8), parameter   :: CON=1.4_8, CON2=1.96_8, BIG=1.0e30_8, SAFE=2.0_8
        integer,      parameter   :: NTAB=10
        integer                   :: i, j
        real(kind=8)              :: errt, fac, hh, a(NTAB,NTAB)
        ! Returns the derivative of a function func at a point x by Ridders’ method of polynomial
        ! extrapolation. The value h is input as an estimated initial stepsize; it need not be small,
        ! but rather should be an increment in x over which func changes substantially. An estimate
        ! of the error in the derivative is returned as err.
        ! Parameters: Stepsize is decreased by CON at each iteration. Max size of tableau is set by
        ! NTAB. Return when error is SAFE worse than the best so far.
        if (h .eq. 0.) then
            dfridr = -10.0_8
        end if

        hh = h
        a(1,1) = (func(x+hh) - func(x-hh)) / (2.0*hh)
        rer = BIG

        do i = 2,NTAB !Successive columns in the Neville tableau will go to smaller stepsizes and higher orders of extrapolation. 
            hh = hh / CON
            a(1,i) = (func(x+hh)-func(x-hh)) / (2.0*hh) !Try new, smaller stepsize.
            fac = CON2
            do j = 2,i !Compute extrapolations of various orders, requiring no new function evaluations.
                a(j,i) = (a(j-1,i)*fac-a(j-1,i-1)) / (fac-1.)
                fac = CON2*fac
                errt = max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
                ! The error strategy is to compare each new extrapolation to one order lower, both at
                ! the present stepsize and the previous one.
                if (errt .le. rer) then !If error is decreased, save the improved answer.
                    rer = errt
                    dfridr = a(j,i)
                end if
            end do

            ! If higher order is worse by a significant factor SAFE, then quit early.
            if (abs(a(i,i)-a(i-1,i-1)) .ge. SAFE * rer) then
                return
            end if 
        end do
        return
    end function dfridr

end function lprimepdf



