function mydtrinv ( func, f, linf, lsup, xinf, xsup, ex, sx, delta, maxit, ier )


! this function is a replacement for dtrinv, too slow...
! f       ex      ex+sx
! target, startg, startd )
!
! linf, lsup, xinf, xsup : UNUSED
!
!******************************************************************************/
!
!  Purpose:
!    mydtrinv computes the quantile of a distribution function func;
!    it uses a simple binary search over the cdf of the cumulative distribution.
!    It should run with fewer function calls than dtrinv.
!
!  Licensing:
!    This code is distributed under the GNU LGPL license. 
!
!  Created:
!    16 March 2020; fortran-version 4 december 2020.
!
!  Author:   Denis Cousineau.
!
!  Parameters:
!    func, f (the quantile, between 0 and 1)
!    linf, lsup, xinf, xsup : unused
!    ex, sx: expected and standard deviation of the distribution
!    delta, maxiter: the precision and maximum number of steps
!    ier: error code (unused)

    implicit none

    !  Function
    real(kind=8) :: mydtrinv

    !  Arguments
    real(kind=8), external    :: func
    real(kind=8), intent(in)  :: f
    logical,      intent(in)  :: linf, lsup
    real(kind=8), intent(in)  :: xinf, xsup, ex, sx, delta
    integer,      intent(in)  :: maxit
    integer,      intent(out) :: ier

    ! Local variables
    real(kind=8)            :: fg, fd, mdl, ptg, ptd
    integer                 :: i, iok
    real(kind=8), parameter :: step=0.5

    ! look for two points on either side of target f
    if (f < func(ex, iok) ) then
        do i = 1, 10
            ptd = ex - (i-1) * step * sx
            ptg = ex - i * step * sx
            if (func(ptg, iok) < f) exit
        end do
    else
        do i = 1, 10
            ptg = ex + (i-1) * step * sx
            ptd = ex + i * step * sx
            if (func(ptd, iok) > f) exit
        end do
    end if

    fg = func(ptg, iok);
    fd = func(ptd, iok);

    do i = 1, maxit
        if (fd - fg < delta ) exit

        mdl = (ptd + ptg)/2
        if (func(mdl, iok) < f) then
            ptg = mdl
        else 
            ptd = mdl
        end if 
        fg = func(ptg, iok)
        fd = func(ptd, iok)
    end do

    mydtrinv = mdl
    
end function mydtrinv
