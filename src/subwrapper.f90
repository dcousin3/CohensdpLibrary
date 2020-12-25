! What follows are wrappers to the fortran functions because 
! R cannot link to functions, only to subroutines
subroutine subklprimesmixturepdf( x, n, obsr, obsd, delta, maxitr, ier, res )
    implicit none
    real(kind=8), intent(in)  :: x, n, obsr, obsd, delta
    integer,      intent(in)  :: maxitr
    integer,      intent(out) :: ier
    real(kind=8), intent(out) :: res
    real(kind=8), external    :: klprimesmixturepdf
    res = klprimesmixturepdf( x, n, obsr, obsd, delta, maxitr, ier )
end subroutine subklprimesmixturepdf

subroutine sublprimecdf( x, n, a, delta, maxitr, ier, res )
    implicit none
    real(kind=8), intent(in)  :: x, n, a, delta
    integer,      intent(in)  :: maxitr
    integer,      intent(out) :: ier
    real(kind=8), intent(out) :: res
    real(kind=8), external    :: lprimecdf
    res = lprimecdf( x, n, a, delta, maxitr, ier )
end subroutine sublprimecdf

subroutine sublprimepdf( x, n, a, delta, maxitr, ier, res )
    implicit none
    real(kind=8), intent(in)  :: x, n, a, delta
    integer,      intent(in)  :: maxitr
    integer,      intent(out) :: ier
    real(kind=8), intent(out) :: res
    real(kind=8), external    :: lprimepdf
    res = lprimepdf( x, n, a, delta, maxitr, ier )
end subroutine sublprimepdf

subroutine sublprimeidf( x, n, a, delta, maxitr, ier, res )
    implicit none
    real(kind=8), intent(in)  :: x, n, a, delta
    integer,      intent(in)  :: maxitr
    integer,      intent(out) :: ier
    real(kind=8), intent(out) :: res
    real(kind=8), external    :: lprimeidf
    res = lprimeidf( x, n, a, delta, maxitr, ier )
end subroutine sublprimeidf

! repeat for kprimepdf and for all cdf, idf