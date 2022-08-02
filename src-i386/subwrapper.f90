! What follows are wrappers to the fortran functions because 
! R cannot link to functions, only to subroutines

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! lsecond
! from Cousineau, 2022, psyArxiv; Lecoutre, 2022
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE sublsecondpdf( delta, n, d, rho, TOL, MAXITER, ier, res )
    IMPLICIT NONE
    INTEGER, PARAMETER         :: PR=KIND(1.0D0)
    REAL(KIND=PR), INTENT(IN)  :: delta, n, d, rho, TOL
    INTEGER,       INTENT(IN)  :: MAXITER
    INTEGER,       INTENT(OUT) :: ier
    REAL(KIND=PR), INTENT(OUT) :: res
    REAL(KIND=PR), EXTERNAL    :: lsecondpdf
    res = lsecondpdf( delta, n, d, rho, TOL, MAXITER, ier )
END SUBROUTINE sublsecondpdf

SUBROUTINE sublsecondcdf( delta, n, d, rho, TOL, MAXITER, ier, res )
    IMPLICIT NONE
    INTEGER, PARAMETER         :: PR=KIND(1.0D0)
    REAL(KIND=PR), INTENT(IN)  :: delta, n, d, rho, TOL
    INTEGER,       INTENT(IN)  :: MAXITER
    INTEGER,       INTENT(OUT) :: ier
    REAL(KIND=PR), INTENT(OUT) :: res
    REAL(KIND=PR), EXTERNAL    :: lsecondcdf
    res = lsecondcdf( delta, n, d, rho, TOL, MAXITER, ier )
END SUBROUTINE sublsecondcdf

SUBROUTINE sublsecondidf( q, n, d, rho, TOL, MAXITER, ier, res )
    IMPLICIT NONE
    INTEGER, PARAMETER         :: PR=KIND(1.0D0)
    REAL(KIND=PR), INTENT(IN)  :: q, n, d, rho, TOL
    INTEGER,       INTENT(IN)  :: MAXITER
    INTEGER,       INTENT(OUT) :: ier
    REAL(KIND=PR), INTENT(OUT) :: res
    REAL(KIND=PR), EXTERNAL    :: lsecondidf
    res = lsecondidf( q, n, d, rho, TOL, MAXITER, ier )
END SUBROUTINE sublsecondidf



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! lprime and kprime distributions
! taken from Poitevineau & Lecoutre, 2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE sublprimepdf( x, nu, a, TOL, MAXITER, ier, res )
    IMPLICIT NONE
    INTEGER, PARAMETER         :: PR=KIND(1.0D0)
    REAL(KIND=PR), INTENT(IN)  :: x, nu, a, TOL
    INTEGER,       INTENT(IN)  :: MAXITER
    INTEGER,       INTENT(OUT) :: ier
    REAL(KIND=PR), INTENT(OUT) :: res
    REAL(KIND=PR), EXTERNAL    :: lprimepdf
    res = lprimepdf( x, nu, a, TOL, MAXITER, ier )
END SUBROUTINE sublprimepdf

SUBROUTINE sublprimecdf( x, nu, a, TOL, MAXITER, ier, res )
    IMPLICIT NONE
    INTEGER, PARAMETER         :: PR=KIND(1.0D0)
    REAL(KIND=PR), INTENT(IN)  :: x, nu, a, TOL
    INTEGER,       INTENT(IN)  :: MAXITER
    INTEGER,       INTENT(OUT) :: ier
    REAL(KIND=PR), INTENT(OUT) :: res
    REAL(KIND=PR), EXTERNAL    :: lprimecdf
    res = lprimecdf( x, nu, a, TOL, MAXITER, ier )
END SUBROUTINE sublprimecdf

SUBROUTINE sublprimeidf( q, nu, a, TOL, MAXITER, ier, res )
    IMPLICIT NONE
    INTEGER, PARAMETER         :: PR=KIND(1.0D0)
    REAL(KIND=PR), INTENT(IN)  :: q, nu, a, TOL
    INTEGER,       INTENT(IN)  :: MAXITER
    INTEGER,       INTENT(OUT) :: ier
    REAL(KIND=PR), INTENT(OUT) :: res
    REAL(KIND=PR), EXTERNAL    :: lprimeidf
    res = lprimeidf( q, nu, a, TOL, MAXITER, ier )
END SUBROUTINE sublprimeidf

SUBROUTINE subkprimepdf( x, nu1, nu2, a, TOL, MAXITER, ier, res )
    IMPLICIT NONE
    INTEGER, PARAMETER         :: PR=KIND(1.0D0)
    REAL(KIND=PR), INTENT(IN)  :: x, nu1, nu2, a, TOL
    INTEGER,       INTENT(IN)  :: MAXITER
    INTEGER,       INTENT(OUT) :: ier
    REAL(KIND=PR), INTENT(OUT) :: res
    REAL(KIND=PR), EXTERNAL    :: kprimepdf
    res = kprimepdf( x, nu1, nu2, a, TOL, MAXITER, ier )
END SUBROUTINE subkprimepdf

SUBROUTINE subkprimecdf( x, nu1, nu2, a, TOL, MAXITER, ier, res )
    IMPLICIT NONE
    INTEGER, PARAMETER         :: PR=KIND(1.0D0)
    REAL(KIND=PR), INTENT(IN)  :: x, nu1, nu2, a, TOL
    INTEGER,       INTENT(IN)  :: MAXITER
    INTEGER,       INTENT(OUT) :: ier
    REAL(KIND=PR), INTENT(OUT) :: res
    REAL(KIND=PR), EXTERNAL    :: kprimecdf
    res = kprimecdf( x, nu1, nu2, a, TOL, MAXITER, ier )
END SUBROUTINE subkprimecdf

SUBROUTINE subkprimeidf( q, nu1, nu2, a, TOL, MAXITER, ier, res )
    IMPLICIT NONE
    INTEGER, PARAMETER         :: PR=KIND(1.0D0)
    REAL(KIND=PR), INTENT(IN)  :: q, nu1, nu2, a, TOL
    INTEGER,       INTENT(IN)  :: MAXITER
    INTEGER,       INTENT(OUT) :: ier
    REAL(KIND=PR), INTENT(OUT) :: res
    REAL(KIND=PR), EXTERNAL    :: kprimeidf
    res = kprimeidf( q, nu1, nu2, a, TOL, MAXITER, ier )
END SUBROUTINE subkprimeidf



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! hypergeometric functions
! from various sources; consult each file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE subhyg0f1( parama, paramz, res )
    IMPLICIT NONE
    INTEGER, PARAMETER         :: PR=KIND(1.0D0)
    REAL(KIND=PR), INTENT(IN)  :: parama, paramz
    REAL(KIND=PR), INTENT(OUT) :: res
    REAL(KIND=PR), EXTERNAL    :: hyg0f1
    res = hyg0f1( parama, paramz )
END SUBROUTINE subhyg0f1

SUBROUTINE subhyg1f1( parama, paramb, paramz, res )
    IMPLICIT NONE
    INTEGER, PARAMETER         :: PR=KIND(1.0D0)
    REAL(KIND=PR), INTENT(IN)  :: parama, paramb, paramz
    REAL(KIND=PR), INTENT(OUT) :: res
    REAL(KIND=PR), EXTERNAL    :: hyg1f1
    res = hyg1f1( parama, paramb, paramz )
END SUBROUTINE subhyg1f1

SUBROUTINE subhyg2f1( parama, paramb, paramc, paramz, res )
    IMPLICIT NONE
    INTEGER, PARAMETER         :: PR=KIND(1.0D0)
    REAL(KIND=PR), INTENT(IN)  :: parama, paramb, paramc, paramz
    REAL(KIND=PR), INTENT(OUT) :: res
    REAL(KIND=PR), EXTERNAL    :: hyg2f1
    res = hyg2f1( parama, paramb, paramc, paramz )
END SUBROUTINE subhyg2f1

