function ncdf( x )

!-----------------------------------------------------------------------
!
!     Returns the probability that a random variable distributed
!     according to the Normal distribution with zero mean and unit
!     variance, is less than or equal to x.
!
!     X    - Input . Argument                                 - Real
!
!     Fortran  functions called:
!        ABS
!
!     Calculation is based upon the error function, using Chebyshev
!     approximation over the interval (0, 6.09).
!
!-----------------------------------------------------------------------

   implicit none

   !  Function
   !  --------

   real(kind=8) :: ncdf

   !  Arguments
   !  ---------

   real(kind=8), intent(in) :: x

   !  Local declarations
   !  ------------------

   real(kind=8), parameter :: zero=0.0_8, half=0.5_8, one=1.0_8
   real(kind=8), parameter :: b=6.09_8
   real(kind=8), parameter :: epsl=1.0e-15_8
   real(kind=8), parameter :: sq2=0.70710678118654748_8  ! = 1/sqrt(2)
   integer, parameter :: m=43

   real(kind=8) :: ax, d, dd, sv, y, y2, z
   integer :: j
   !  Chebyshev coefficients for the error function:
   real(kind=8) :: c(m)=(/                                          &
                0.1635454272828649e+01_8, 0.3340345052068802e+00_8, &
               -0.2550158252490571e+00_8, 0.1576322031081386e+00_8, &
               -0.7292374105725899e-01_8, 0.1844981586051652e-01_8, &
                0.5451830906274930e-02_8,-0.9378172345987997e-02_8, &
                0.5461774234132427e-02_8,-0.1369434835987556e-02_8, &
               -0.4791673249620534e-03_8, 0.6448576838517482e-03_8, &
               -0.2787148092871001e-03_8, 0.1375372712760661e-04_8, &
                0.5421909937317043e-04_8,-0.3178226588449582e-04_8, &
                0.4982720716690025e-05_8, 0.3910358552338792e-05_8, &
               -0.2794516331451236e-05_8, 0.5298695171317318e-06_8, &
                0.2745633118680364e-06_8,-0.2071635016913643e-06_8, &
                0.3687897944251272e-07_8, 0.1967578278235985e-07_8, &
               -0.1324655993569593e-07_8, 0.1759666490397310e-08_8, &
                0.1371495927659758e-08_8,-0.7222588259009134e-09_8, &
                0.4293513190125905e-10_8, 0.8604672660411214e-10_8, &
               -0.3240778980162313e-10_8,-0.1531203125775596e-11_8, &
                0.4585991808162509e-11_8,-0.1112803317375455e-11_8, &
               -0.2501311672927168e-12_8, 0.1996065316425628e-12_8, &
               -0.2333998072735281e-13_8,-0.1670821608490672e-13_8, &
                0.6679068635371678e-14_8, 0.2100528246172614e-15_8, &
               -0.7805650505799361e-15_8,-0.2386236005000803e-15_8, &
               -0.1792153097738274e-15_8 /)

   !--------------------------------------------------------------------

   if ( abs(x) < epsl ) then
      ncdf = half
      return
   end if

   z  = x * sq2
   ax = abs( z )
   if ( ax < b ) then
      d  = zero
      dd = zero
      y  = (ax+ax-b)/b
      y2 = y + y
      do j = m, 2, -1
         sv = d
         d  = y2*d - dd + c(j)
         dd = sv
      end do
      ncdf = half + half*(y*d-dd+half*c(1))
   else
      ncdf = one
   end if
   if ( x < zero ) ncdf = one - ncdf

end function ncdf
