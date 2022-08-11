FUNCTION FBlsecondCDF( delta, n, d, r, rho, TOL, MAXITER )

    !! This distribution is based on a fiducial-Bayesian approach
    !! in which the mixing parameter K is based on an observed r rather than the population rho.

    ! ***************************************************************************
    ! Nonstandard noncentral Lambda" distribution
    !   This distribution is the exact distribution of the Cohen's delta
    !   in the 2-repeated-measure samples from a Binormal population with 
    !   common variance sigma.
    ! Parameters are n        (sample size), 
    !                d        (observed difference in the sample), 
    !                r        (observed correlation in the sample)
    !                rho      (correlation in the population)
    !                delta    (putative parameter whose probability is required)
    !                MAXITER  (maximum iterations, often 15 are enough)
    !                TOL      (a tolerance level to decide to quit)
    ! Outputs the probability of delta given the population parameter rho and observed d_p
    !
    ! Requires lprimecdf, LOG_GAMMA, EXP
    !
    ! Because there are many large numbers, this is far from trivial to have
    !    a robust program...
    !
    ! Denis Cousineau 12/juillet/2022.
    !
    ! ***************************************************************************
    IMPLICIT NONE
    INTEGER, PARAMETER      :: PR=KIND(1.0D0)

    REAL(PR), intent(in)    :: delta, n, d, r, rho
    integer                 :: MAXITER  ! max iterations, suggested 5000
    REAL(PR), INTENT(IN)    :: TOL      ! precision to exit, suggested 10D-5
    REAL(PR)                :: FBlsecondCDF

    REAL(PR), EXTERNAL      :: lprimecdf
    REAL(PR), PARAMETER     :: HALF=0.5D0,ONE=1.0D0,TWO=2.0D0,FOUR=4.0D0,SEIZE=16.0D0

    !  Local declarations
    REAL(PR)              :: q, b0, total, step, previousstep, ck, ek, ier
    INTEGER               :: k

    step  = 0.0D0
    previousstep = 0.0D0
    total  = 0.0D0

    b0 = sqrt(n/(TWO*(ONE-rho))) 
    q  = n - ONE

    k = 0
    DO WHILE ( k .LE. MAXITER )
        previousstep = step
        ! deux changements ici q -> q-1 (3x) et rho -> r (2x)
        ck =  EXP( LOG_GAMMA((q-1)/TWO+k) - LOG_GAMMA((q-1)/TWO) - LOG_GAMMA(k+ONE) + &
                k * LOG(r**2) + ((q-1)/TWO) * LOG(ONE-r**2)                           &
        )
        if (r .EQ. 0.0D0) ck = 0.0D0
        if (r .EQ. 1.0D0) ck = 0.0D0

        ek = sqrt(1/(1-rho**2)) * sqrt(q/(q + TWO*k))

        step = ck * lprimecdf(b0 * delta, TWO*q+FOUR*k, b0*d / ek, TOL, MAXITER, ier)

        total = total +  step
        k = k + 1
        ! exit if too many iterations or if step is decreasing and goes below TOL
        IF (( k > MAXITER).or.((step < previousstep).AND.(step < TOL))) GOTO 10

    ENDDO

10  FBlsecondCDF =  total

END FUNCTION FBlsecondCDF

