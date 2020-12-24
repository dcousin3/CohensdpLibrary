function klprimesmixturecdf( x, n, obsr, obsd, delta, maxitr, ier )

!-----------------------------------------------------------------------
!
!     Calculates the probability that a random variable distributed
!     according to the L' distribution with parameters 
!           n        degree of freedom, and
!           obsd     centrality parameter itself a function of rho 
!     where rho is distributed according to the K' distribution with parameters
!           n-2, n-1 degrees of freedom and
!           obsr     centrality parameter
!     is less than or equal to X. This distribution is a convolution across all rhos
!     given one observed r. It is expressed as
!           CDF(x) = NIntegrate[PDF[rDistribution[n-2, n-1, obsr], rho] *
!                               CDF[dDistribution[n, obsd, rho], x], 
!                               {rho, -1, 1}]
!     where
!            PDF[rDistribution[n-2, n-1, obsr], rho] = abs(t') PDF[K'(n-2,n-1,sqrt(n-2] rho/sqrt(1-rho^2)](sqrt(n-1) r/sqrt(1-r^2))
!            CDF[dDistribution[n,dp, rho], delta] = sqrt(n/(2(1-rho)) PDF[L'(2(n-1)/(1+rho^2), sqrt(n/(2(1-rho))) dp](sqrt(n/(2(1-rho))) x)
!
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
!       KPRIMECDF  LPRIMECDF
!     Fortran functions called:
!       ABS  LOG  SQRT  DLGAMA
!
!*********************************************************************************************!
!**                                                                                         **!
!** This function was created by Denis Cousineau, 3 december 2020.                          **!
!**                                                                                         **!
!** It uses functions from B. Lecoutre, programmed by J. Poitevineau.                       **!
!** It also uses functions from  Press, Teukolsky, Vetterling and Flanery, 1992             **!
!**                                                                                         **!
!*********************************************************************************************!

    implicit none

    !  Function
    real(kind=8) :: klprimesmixturecdf

    !  Arguments
    real(kind=8), intent(in)  :: x, n, obsr, obsd, delta
    integer,      intent(in)  :: maxitr
    integer,      intent(out) :: ier
    real(kind=8), parameter        :: EPS=1.e-6_8

    !  Local declarations
    real(kind=8)              :: s
    real(kind=8), external    :: kprimepdf, lprimecdf

    ! NOTE: the bounds cannot be exactly -1 and +1 so removed EPS
    s = qromb( fct, -1.0_8 + EPS, +1.0_8 - EPS )

    ! s is the solution
    ier = 0
    klprimesmixturecdf = s

contains
    function fct(rho)
        real(kind=8), intent(in)  :: rho
        real(kind=8), external    :: kprimepdf, lprimecdf
        real(kind=8)              :: fct

        real(kind=8)              :: a, dfL, ncL, ncK, scaledr, scaledd !temporary variables
        integer                   :: iok

        a = abs( (sqrt(n-1._8) * rho**2 )/((1-rho**2)**(1.5_8)) + sqrt(n-1._8)/sqrt(1-rho**2) )
        scaledr = sqrt(n-1._8) * rho  / sqrt(1-rho**2)
        ncK     = sqrt(n-2._8) * obsr / sqrt(1-obsr**2)

        dfL     = 2 * (n-1) / (1+rho**2)             ! Cousineau, 2020
        ncL     = sqrt(n / (2 * (1 - rho) ) ) * obsd ! Becker, 1988
        scaledd = sqrt(n / (2 * (1 - rho) ) ) * x    ! Becker, 1988

        fct = a  * kprimepdf(scaledr, n-2.0_8, n-1.0_8, ncK, delta, maxitr, iok)  &
                 * lprimecdf(scaledd, dfL, ncL, delta, maxitr, iok)

    end function fct

    ! The following functions are taken from Press, Teukolsky, Vetterling, and Flannery
    ! Numerical Receipes in Fortran 77, 1992.
    function qromb(fct, a, b)
        real(kind=8), external         :: fct
        real(kind=8), intent(in)       :: a, b
        real(kind=8)                   :: qromb

        real(kind=8), parameter        :: EPS=1.e-6_8
        integer, parameter             :: JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1
        integer                        :: j
        real(kind=8)                   :: dqromb
        real(kind=8), dimension(JMAXP) :: h, s 

        h(1) = 1.0_8
        do j = 1, JMAX
            call trapzd(fct, a, b, s(j), j)
            if (j .ge. K) then
                call polint(h(j-KM), s(j-KM), K, 0.0_8, qromb, dqromb)
                if (abs(dqromb) .le. EPS*abs(qromb)) return
            endif
            s(j+1) = s(j)
            h(j+1) = 0.25*h(j)
        enddo
        qromb = -43.0_8 ! something went wrong
        return
    end function qromb


    subroutine trapzd(fct, a, b, s, n)
        integer, intent(in)         :: n
        real(kind=8), external      :: fct
        real(kind=8), intent(in)    :: a, b
        real(kind=8), intent(inout) :: s

        integer                     :: it, j
        real(kind=8)                :: del, fsum, tnm, x

        if (n .eq. 1) then
            s = 0.5*(b-a)*(fct(a)+fct(b))
        else
            it  = 2**(n-2)
            tnm = it
            del = (b-a)/tnm 
            x   = a+0.5*del
            fsum = 0.
            do j = 1, it
                fsum = fsum + fct(x)
                x   = x + del
            enddo
            s = 0.5*(s+(b-a)*fsum/tnm) 
        endif
        return
    end subroutine trapzd


    subroutine polint(xa, ya, n, x, y, dy)
        integer, intent(in)        :: n
        real(kind=8), intent(in)   :: xa(n), ya(n), x
        real(kind=8), intent(out)  :: y, dy

        integer, parameter         :: NMAX=10
        integer                    :: i, m, ns
        real(kind=8)               :: den, dif, dift, ho, hp, w, c(NMAX), d(NMAX)

        ns = 1
        dif = abs(x-xa(1))
        do i = 1, n 
            dift = abs(x-xa(i))
            if (dift .lt. dif) then
                ns = i
                dif = dift
            endif
            c(i) = ya(i) 
            d(i) = ya(i)
        enddo
        y = ya(ns) 
        ns = ns-1
        do m = 1, n-1 
            do i = 1, n-m 
                ho   = xa(i)-x
                hp   = xa(i+m)-x
                w    = c(i+1)-d(i)
                den  = ho-hp
                if (den .eq. 0.) then
                    y  = -45.0_8
                    dy = -45.0_8
                    return ! pause "failure in polint"
                endif
                den  = w/den
                d(i) = hp*den 
                c(i) = ho*den
            enddo
            if (2*ns .lt. n-m) then 
                dy = c(ns+1)
            else
                dy = d(ns)
                ns = ns - 1
            endif
            y = y + dy
        enddo
        return
    end subroutine polint


end function klprimesmixturecdf

