!      ALGORITHM 707, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 3, SEPTEMBER, 1992, PP. 345-349.
!     ****************************************************************  
!     *      SOLUTION TO THE CONFLUENT HYPERGEOMETRIC FUNCTION       *  
!     *                           by                                 *  
!     *                      MARK NARDIN,                            *  
!     *              W. F. PERGER and ATUL BHALLA                    *  
!     *  Michigan Technological University, Copyright 1989           *  
!     *  Adapted 04/04/2022 Denis Cousineau for F90                  *  
!     *                                                              *  
!     *  Description : A numerical evaluator for the confluent       *  
!     *    hypergeometric function for complex arguments with large  *  
!     *    magnitudes using a direct summation of the Kummer series. *  
!     *    The method used allows an accuracy of up to thirteen      *  
!     *    decimal places through the use of large real arrays       *  
!     *    and a single final division.  LNCHF is a variable which   *  
!     *    selects how the result should be represented.  A '0' will *  
!     *    return the value in standard exponential form.  A '1'     *  
!     *    will return the LOG of the result.  IP is an integer      *
!     *    variable that specifies how many array positions are      *
!     *    desired (usually 10 is sufficient).  Setting IP=0 causes  *
!     *    the program to estimate the number of array positions.    *
!     *                                                              *  
!     *    The confluent hypergeometric function is the solution to  *  
!     *    the differential equation:                                *  
!     *             zf"(z) + (a-z)f'(z) - bf(z) = 0                  *  
!     *  Subprograms called: BITS, CHGF                              *  
!     ****************************************************************  

! addition for real numbers only by D. Cousineau, 4/4/2022 (at the begining).

FUNCTION hyg1F1(A, B, Z)
    INTEGER, PARAMETER    :: PR=KIND(1.0D0)
    REAL(PR), INTENT(IN)  :: A, B, Z
    REAL(PR)              :: hyg1F1
!    COMPLEX(PR), EXTERNAL :: CONHYP
    REAL(PR)              :: HG

    ! unused: precision is too poor...
    ! hyg1F1 = REAL( CONHYP(CMPLX(A, 0, KIND=PR),CMPLX(B, 0, KIND=PR), & 
    !                     CMPLX(Z, 0, KIND=PR), 0, 200 ) )

    ! Moreau version, more precise.
    CALL CHGM(A,B,Z,HG)
    hyg1F1 = HG

END FUNCTION hyg1F1



    !******************************************************************
    !*      Purpose: This program computes the confluent              *
    !*               hypergeometric function M(a,b,x) using           *
    !*               subroutine CHGM                                  *
    !*      Input  : a  --- Parameter                                 *
    !*               b  --- Parameter ( b <> 0,-1,-2,... )            *
    !*               x  --- Argument                                  *
    !*      Output:  HG --- M(a,b,x)                                  *
    !* -------------------------------------------------------------- *
    !* REFERENCE: "Fortran Routines for Computation of Special        *
    !*             Functions jin.ece.uiuc.edu/routines/routines.html" *
    !*                                                                *
    !*                              F90 Release By J-P Moreau, Paris. *
    !*                                     (www.jpmoreau.fr)          *
    !******************************************************************

SUBROUTINE CHGM(A,B,X,HG)
    !===================================================
    !  Purpose: Compute confluent hypergeometric function
    !           M(a,b,x)
    !  Input  : a  --- Parameter
    !           b  --- Parameter ( b <> 0,-1,-2,... )
    !           x  --- Argument
    !  Output:  HG --- M(a,b,x)
    !  Routine called: GAMMA for computing ?(x)
    !===================================================
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    PI=3.141592653589793D0
    A0=A
    A1=A
    X0=X
    HG=0.0D0
    IF (B.EQ.0.0D0.OR.B.EQ.-ABS(INT(B))) THEN
       HG=1.0D+300
    ELSE IF (A.EQ.0.0D0.OR.X.EQ.0.0D0) THEN
       HG=1.0D0
    ELSE IF (A.EQ.-1.0D0) THEN
       HG=1.0D0-X/B
    ELSE IF (A.EQ.B) THEN
       HG=DEXP(X)
    ELSE IF (A-B.EQ.1.0D0) THEN
       HG=(1.0D0+X/B)*DEXP(X)
    ELSE IF (A.EQ.1.0D0.AND.B.EQ.2.0D0) THEN
       HG=(DEXP(X)-1.0D0)/X
    ELSE IF (A.EQ.INT(A).AND.A.LT.0.0D0) THEN
       M=INT(-A)
       R=1.0D0
       HG=1.0D0
       DO 10 K=1,M
          R=R*(A+K-1.0D0)/K/(B+K-1.0D0)*X
10        HG=HG+R
    ENDIF
    IF (HG.NE.0.0D0) RETURN
    IF (X.LT.0.0D0) THEN
       A=B-A
       A0=A
       X=DABS(X)
    ENDIF
    IF (A.LT.2.0D0) NL=0
    IF (A.GE.2.0D0) THEN
       NL=1
       LA=INT(A)
       A=A-LA-1.0D0
    ENDIF
    DO 30 N=0,NL
       IF (A0.GE.2.0D0) A=A+1.0D0
       IF (X.LE.30.0D0+DABS(B).OR.A.LT.0.0D0) THEN
          HG=1.0D0
          RG=1.0D0
          DO 15 J=1,500
             RG=RG*(A+J-1.0D0)/(J*(B+J-1.0D0))*X
             HG=HG+RG
             IF (DABS(RG/HG).LT.1.0D-15) GO TO 25
15        CONTINUE
       ELSE
          TA =  GAMMA(A)
          TB =  GAMMA(B)
          XG=B-A
          TBA = GAMMA(XG)
          SUM1=1.0D0
          SUM2=1.0D0
          R1=1.0D0
          R2=1.0D0
          DO 20 I=1,8
             R1=-R1*(A+I-1.0D0)*(A-B+I)/(X*I)
             R2=-R2*(B-A+I-1.0D0)*(A-I)/(X*I)
             SUM1=SUM1+R1
20               SUM2=SUM2+R2
          HG1=TB/TBA*X**(-A)*DCOS(PI*A)*SUM1
          HG2=TB/TA*DEXP(X)*X**(A-B)*SUM2
          HG=HG1+HG2
       ENDIF
25         IF (N.EQ.0) Y0=HG
       IF (N.EQ.1) Y1=HG
30      CONTINUE
    IF (A0.GE.2.0D0) THEN
       DO 35 I=1,LA-1
          HG=((2.0D0*A-B+X)*Y1+(B-A)*Y0)/A
          Y0=Y1
          Y1=HG
35        A=A+1.0D0
    ENDIF

    IF (X0.LT.0.0D0) HG=HG*DEXP(X0)
    A=A1
    X=X0
    RETURN
END






!! All that follow not USED: precision is too poor...


FUNCTION CONHYP (A,B,Z,LNCHF,IP)                          
                                                                        
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)
      INTEGER         :: LNCHF,I,BITS,IP
      COMPLEX(KIND=PR) :: CHGF,A,B,Z,CONHYP                                      
      REAL(KIND=PR)    :: NTERM,FX,TERM1,MAX,TERM2,ANG                               
                                                                        
      IF (ABS(Z) .NE. 0.0D0) THEN                                     
        ANG=ATAN2(AIMAG(Z),REAL(Z, KIND=PR))                                   
      ELSE                                                              
        ANG=1.0D0                                                       
      ENDIF                                                             
      IF (DABS(ANG) .LT. (3.14159D0*0.5)) THEN                          
        ANG=1.0D0                                                       
      ELSE                                                              
        ANG=DSIN(DABS(ANG)-(3.14159265D0*0.5D0))+1.0D0                  
      ENDIF                                                             
      MAX=0                                                             
      NTERM=0                                                           
      FX=0                                                              
      TERM1=0                                                           
10    NTERM=NTERM+1                                                     
      TERM2=ABS((A+NTERM-1)*Z/((B+NTERM-1)*NTERM))                    
      IF (TERM2 .EQ. 0.0D0) GOTO 20                                     
      IF (TERM2 .LT. 1.0D0) THEN                                        
        IF ((DBLE(A)+NTERM-1) .GT. 1.0D0) THEN                         
          IF ((DBLE(B)+NTERM-1) .GT. 1.0D0) THEN                       
            IF ((TERM2-TERM1) .LT. 0.0D0) THEN                          
              GOTO 20                                                   
            ENDIF                                                       
          ENDIF                                                         
        ENDIF                                                           
      ENDIF                                                             
      FX=FX+DLOG(TERM2)                                                 
      IF (FX .GT. MAX) MAX=FX                                           
      TERM1=TERM2                                                       
      GOTO 10                                                           
20    MAX=MAX*2/(BITS()*6.93147181D-1)                                  
      I=INT(MAX*ANG)+7                                                  
      IF (I .LT. 5) I=5                                                 
      IF (IP .GT. I) I=IP
      CONHYP=CHGF(A,B,Z,I,LNCHF)                                        
      RETURN                                                            
END FUNCTION CONHYP                                                             


                                                                        
INTEGER FUNCTION BITS()                                           
!     ****************************************************************
!     *                   FUNCTION BITS                              *  
!     *                                                              *  
!     *                                                              *  
!     *  Description : Determines the number of significant figures  *  
!     *    of machine precision to arrive at the size of the array   *  
!     *    the numbers must must be stored in to get the accuracy    *  
!     *    of the solution.                                          *  
!     *  Subprogram called: STORE                                    *
!     ****************************************************************  
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)
      REAL(KIND=PR) :: BIT,BIT2,STORE
      INTEGER      :: COUNT                                                     
                                                                        
      BIT   = 1.0                                                           
      COUNT = 0                                                           
10    COUNT = COUNT+1                                                     
      BIT2  = STORE(BIT*2.0)
      BIT   = STORE(BIT2+1.0)
      IF ((BIT-BIT2) .NE. 0.0) GOTO 10                                  
      BITS=COUNT
      RETURN                                                            
END FUNCTION BITS


DOUBLE PRECISION FUNCTION STORE (X)
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)
      REAL(KIND=PR) ::  X
!***********************************************************
!   This function forces its argument X to be stored in a
! memory location, thus providing a means of determining
! floating point number characteristics (such as the machine
! precision) when it is necessary to avoid computation in
! high precision registers.
! On input:
!       X = Value to be stored.
!       X is not altered by this function.
! On output:
!       STORE = Value of X after it has been stored and
!               possibly truncated or rounded to the double
!               precision word length.
! Modules required by STORE:  None
!***********************************************************
      REAL(KIND=PR) :: Y
      COMMON/STCOM/Y
      Y = X
      STORE = Y
      RETURN
END FUNCTION STORE
                                                                        
                                                                        
FUNCTION CHGF (A,B,Z,L,LNCHF)                          
!     ****************************************************************  
!     *                   FUNCTION CHGF                              *  
!     *  Description : Function that sums the Kummer series and      *  
!     *    returns the solution of the confluent hypergeometric      *  
!     *    function.                                                 *  
!     *  Subprograms called: ARMULT, ARYDIV, BITS, CMPADD, CMPMUL    *  
!     ****************************************************************  
      PARAMETER (LENGTH=25000)
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)
      INTEGER         :: L,I,BITS,BIT,LNCHF
      COMPLEX(KIND=PR) :: A,B,Z,FINAL,CHGF                                       
      REAL(KIND=PR) :: AR,AI,CR,CI,XR,XI,CNT,SIGFIG,MX1,MX2,RMAX
      REAL(KIND=PR), DIMENSION(-1:LENGTH) :: SUMR,SUMI,NUMR,NUMI,DENOMR,DENOMI, QR1,QR2,QI1,QI2
      REAL(KIND=PR) :: AR2,AI2,CR2,CI2,XR2,XI2
 
                                                                        
      BIT=BITS()                                                        
      RMAX=2.0D0**(BIT/2)                                               
      SIGFIG=2.0D0**(BIT/4)                                             
      AR2=DBLE(A)*SIGFIG                                               
      AR=DINT(AR2)                                                      
      AR2=DNINT((AR2-AR)*RMAX)                                          
      AI2=AIMAG(A)*SIGFIG                                               
      AI=DINT(AI2)                                                      
      AI2=DNINT((AI2-AI)*RMAX)                                          
      CR2=DBLE(B)*SIGFIG                                               
      CR=DINT(CR2)                                                      
      CR2=DNINT((CR2-CR)*RMAX)                                          
      CI2=AIMAG(B)*SIGFIG                                               
      CI=DINT(CI2)                                                      
      CI2=DNINT((CI2-CI)*RMAX)                                          
      XR2=DBLE(Z)*SIGFIG                                               
      XR=DINT(XR2)                                                      
      XR2=DNINT((XR2-XR)*RMAX)                                          
      XI2=AIMAG(Z)*SIGFIG                                               
      XI=DINT(XI2)                                                      
      XI2=DNINT((XI2-XI)*RMAX)                                          
      SUMR(-1)=1.0D0                                                    
      SUMI(-1)=1.0D0                                                    
      NUMR(-1)=1.0D0                                                    
      NUMI(-1)=1.0D0                                                    
      DENOMR(-1)=1.0D0                                                  
      DENOMI(-1)=1.0D0                                                  
      DO 100 I=0,L+1                                                    
        SUMR(I)=0.0D0                                                   
        SUMI(I)=0.0D0                                                   
        NUMR(I)=0.0D0                                                   
        NUMI(I)=0.0D0                                                   
        DENOMR(I)=0.0D0                                                 
        DENOMI(I)=0.0D0                                                 
100   CONTINUE                                                          
      SUMR(1)=1.0D0                                                     
      NUMR(1)=1.0D0                                                     
      DENOMR(1)=1.0D0                                                   
      CNT=SIGFIG                                                        
110   IF (SUMR(1) .LT. 0.5) THEN                                        
        MX1=SUMI(L+1)                                                   
      ELSE IF (SUMI(1) .LT. 0.5) THEN                                   
        MX1=SUMR(L+1)                                                   
      ELSE                                                              
        MX1=DMAX1(SUMR(L+1),SUMI(L+1))                                  
      ENDIF                                                             
      IF (NUMR(1) .LT. 0.5) THEN                                        
        MX2=NUMI(L+1)                                                   
      ELSE IF (NUMI(1) .LT. 0.5) THEN                                   
        MX2=NUMR(L+1)                                                   
      ELSE                                                              
        MX2=DMAX1(NUMR(L+1),NUMI(L+1))                                  
      ENDIF                                                             
      IF (MX1-MX2 .GT.  2.0) THEN                                       
        IF (CR .GT. 0.0D0) THEN                                         
          IF (ABS(CMPLX(AR,AI)*CMPLX(XR,XI)/(CMPLX(CR,CI)*CNT)) .LE. 1.0D0) GOTO 190
        ENDIF                                                           
      ENDIF                                                             
      CALL CMPMUL(SUMR,SUMI,CR,CI,QR1,QI1,L,RMAX)                       
      CALL CMPMUL(SUMR,SUMI,CR2,CI2,QR2,QI2,L,RMAX)                     
      QR2(L+1)=QR2(L+1)-1                                               
      QI2(L+1)=QI2(L+1)-1                                               
      CALL CMPADD(QR1,QI1,QR2,QI2,SUMR,SUMI,L,RMAX)                     
                                                                        
      CALL ARMULT(SUMR,CNT,SUMR,L,RMAX)                                 
      CALL ARMULT(SUMI,CNT,SUMI,L,RMAX)                                 
      CALL CMPMUL(DENOMR,DENOMI,CR,CI,QR1,QI1,L,RMAX)                   
      CALL CMPMUL(DENOMR,DENOMI,CR2,CI2,QR2,QI2,L,RMAX)                 
      QR2(L+1)=QR2(L+1)-1                                               
      QI2(L+1)=QI2(L+1)-1                                               
      CALL CMPADD(QR1,QI1,QR2,QI2,DENOMR,DENOMI,L,RMAX)                 
                                                                        
      CALL ARMULT(DENOMR,CNT,DENOMR,L,RMAX)                             
      CALL ARMULT(DENOMI,CNT,DENOMI,L,RMAX)                             
      CALL CMPMUL(NUMR,NUMI,AR,AI,QR1,QI1,L,RMAX)                       
      CALL CMPMUL(NUMR,NUMI,AR2,AI2,QR2,QI2,L,RMAX)                     
      QR2(L+1)=QR2(L+1)-1                                               
      QI2(L+1)=QI2(L+1)-1                                               
      CALL CMPADD(QR1,QI1,QR2,QI2,NUMR,NUMI,L,RMAX)                     
                                                                        
      CALL CMPMUL(NUMR,NUMI,XR,XI,QR1,QI1,L,RMAX)                       
      CALL CMPMUL(NUMR,NUMI,XR2,XI2,QR2,QI2,L,RMAX)                     
      QR2(L+1)=QR2(L+1)-1                                               
      QI2(L+1)=QI2(L+1)-1                                               
      CALL CMPADD(QR1,QI1,QR2,QI2,NUMR,NUMI,L,RMAX)                     
                                                                        
      CALL CMPADD(SUMR,SUMI,NUMR,NUMI,SUMR,SUMI,L,RMAX)                 
      CNT=CNT+SIGFIG                                                    
      AR=AR+SIGFIG                                                      
      CR=CR+SIGFIG                                                      
      GOTO 110                                                          
190   CALL ARYDIV(SUMR,SUMI,DENOMR,DENOMI,FINAL,L,LNCHF,RMAX,BIT)       
      CHGF=FINAL                                                        
      RETURN                                                            
END FUNCTION CHGF
                                                                        
                                                                        
                                                                        
SUBROUTINE ARADD(A,B,C,L,RMAX)
!     ****************************************************************  
!     *                 SUBROUTINE ARADD                             *  
!     *  Description : Accepts two arrays of numbers and returns     *  
!     *    the sum of the array.  Each array is holding the value    *  
!     *    of one number in the series.  The parameter L is the      *  
!     *    size of the array representing the number and RMAX is     *  
!     *    the actual number of digits needed to give the numbers    *  
!     *    the desired accuracy.                                     *  
!     *  Subprograms called: none                                    *  
!     ****************************************************************  
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)
      INTEGER      :: L
      REAL(KIND=PR) :: A,B,C,Z,RMAX                                               
      INTEGER      :: EDIFF,I,J                                                 

      DIMENSION A(-1:*),B(-1:*),C(-1:*),Z(-1:25000)                 
                                                                        
      DO 110 I=0,L+1                                                    
        Z(I)=0.0D0                                                      
110   CONTINUE                                                          
      EDIFF=DNINT(A(L+1)-B(L+1))                                        
      IF (DABS(A(1)) .LT. 0.5 .OR. EDIFF .LE. -L) GOTO 111              
      IF (DABS(B(1)) .LT. 0.5 .OR. EDIFF .GE. L) GOTO 113               
      GOTO 115                                                          
111   DO 112 I=-1,L+1                                                   
        C(I)=B(I)                                                       
112   CONTINUE                                                          
      GOTO 311                                                          
113   DO 114 I=-1,L+1                                                   
        C(I)=A(I)                                                       
114   CONTINUE                                                          
      GOTO 311                                                          
115   Z(-1)=A(-1)                                                       
      IF (DABS(A(-1)-B(-1)) .LT. 0.5) GOTO 200                          
      IF (EDIFF .GT. 0) THEN                                            
        Z(L+1)=A(L+1)                                                   
        GOTO 233                                                        
      ENDIF                                                             
      IF (EDIFF .LT. 0) THEN                                            
        Z(L+1)=B(L+1)                                                   
        Z(-1)=B(-1)                                                     
        GOTO 266                                                        
      ENDIF                                                             
      DO 120 I=1,L                                                      
        IF (A(I) .GT. B(I)) THEN                                        
          Z(L+1)=A(L+1)                                                 
          GOTO 233                                                      
        ENDIF                                                           
        IF (A(I) .LT. B(I)) THEN                                        
          Z(L+1)=B(L+1)                                                 
          Z(-1)=B(-1)                                                   
          GOTO 266                                                      
        ENDIF                                                           
120   CONTINUE                                                          
      GOTO 300                                                          
                                                                        
200   IF (EDIFF .GT. 0) GOTO 203                                        
      IF (EDIFF .LT. 0) GOTO 207                                        
      Z(L+1)=A(L+1)                                                     
      DO 201 I=L,1,-1                                                   
        Z(I)=A(I)+B(I)+Z(I)                                             
        IF (Z(I) .GE. RMAX) THEN                                        
          Z(I)=Z(I)-RMAX                                                
          Z(I-1)=1.0D0                                                  
        ENDIF                                                           
201   CONTINUE                                                          
      IF (Z(0) .GT. 0.5) THEN                                           
        DO 202 I=L,1,-1                                                 
          Z(I)=Z(I-1)                                                   
202     CONTINUE                                                        
        Z(L+1)=Z(L+1)+1.0D0                                             
        Z(0)=0.0D0                                                      
      ENDIF                                                             
      GOTO 300                                                          
203   Z(L+1)=A(L+1)                                                     
      DO 204 I=L,1+EDIFF,-1                                             
        Z(I)=A(I)+B(I-EDIFF)+Z(I)                                       
        IF (Z(I) .GE. RMAX) THEN                                        
          Z(I)=Z(I)-RMAX                                                
          Z(I-1)=1.0D0                                                  
        ENDIF                                                           
204   CONTINUE                                                          
      DO 205 I=EDIFF,1,-1                                               
        Z(I)=A(I)+Z(I)                                                  
        IF (Z(I) .GE. RMAX) THEN                                        
          Z(I)=Z(I)-RMAX                                                
          Z(I-1)=1.0D0                                                  
        ENDIF                                                           
205   CONTINUE                                                          
      IF (Z(0) .GT. 0.5) THEN                                           
        DO 206 I=L,1,-1                                                 
          Z(I)=Z(I-1)                                                   
206     CONTINUE                                                        
        Z(L+1)=Z(L+1)+1                                                 
        Z(0)=0.0D0                                                      
      ENDIF                                                             
      GOTO 300                                                          
207   Z(L+1)=B(L+1)                                                     
      DO 208 I=L,1-EDIFF,-1                                             
        Z(I)=A(I+EDIFF)+B(I)+Z(I)                                       
        IF (Z(I) .GE. RMAX) THEN                                        
          Z(I)=Z(I)-RMAX                                                
          Z(I-1)=1.0D0                                                  
        ENDIF                                                           
208   CONTINUE                                                          
      DO 209 I=0-EDIFF,1,-1                                             
        Z(I)=B(I)+Z(I)                                                  
        IF (Z(I) .GE. RMAX) THEN                                        
          Z(I)=Z(I)-RMAX                                                
          Z(I-1)=1.0D0                                                  
        ENDIF                                                           
209   CONTINUE                                                          
      IF (Z(0) .GT. 0.5) THEN                                           
        DO 210 I=L,1,-1                                                 
          Z(I)=Z(I-1)                                                   
210     CONTINUE                                                        
        Z(L+1)=Z(L+1)+1.0D0                                             
        Z(0)=0.0D0                                                      
      ENDIF                                                             
      GOTO 300                                                          
                                                                        
233   IF (EDIFF .GT. 0) GOTO 243                                        
      DO 234 I=L,1,-1                                                   
        Z(I)=A(I)-B(I)+Z(I)                                             
        IF (Z(I) .LT. 0.0D0) THEN                                       
          Z(I)=Z(I)+RMAX                                                
          Z(I-1)=-1.0D0                                                 
        ENDIF                                                           
234   CONTINUE                                                          
      GOTO 290                                                          
243   DO 244 I=L,1+EDIFF,-1                                             
        Z(I)=A(I)-B(I-EDIFF)+Z(I)                                       
        IF (Z(I) .LT. 0.0D0) THEN                                       
          Z(I)=Z(I)+RMAX                                                
          Z(I-1)=-1.0D0                                                 
        ENDIF                                                           
244   CONTINUE                                                          
      DO 245 I=EDIFF,1,-1                                               
        Z(I)=A(I)+Z(I)                                                  
        IF (Z(I) .LT. 0.0D0) THEN                                       
          Z(I)=Z(I)+RMAX                                                
          Z(I-1)=-1.0D0                                                 
        ENDIF                                                           
245   CONTINUE                                                          
      GOTO 290                                                          
                                                                        
266   IF (EDIFF .LT. 0) GOTO 276                                        
      DO 267 I=L,1,-1                                                   
        Z(I)=B(I)-A(I)+Z(I)                                             
        IF (Z(I) .LT. 0.0D0) THEN                                       
          Z(I)=Z(I)+RMAX                                                
          Z(I-1)=-1.0D0                                                 
        ENDIF                                                           
267   CONTINUE                                                          
      GOTO 290                                                          
276   DO 277 I=L,1-EDIFF,-1                                             
        Z(I)=B(I)-A(I+EDIFF)+Z(I)                                       
        IF (Z(I) .LT. 0.0D0) THEN                                       
          Z(I)=Z(I)+RMAX                                                
          Z(I-1)=-1.0D0                                                 
        ENDIF                                                           
277   CONTINUE                                                          
      DO 278 I=0-EDIFF,1,-1                                             
        Z(I)=B(I)+Z(I)                                                  
        IF (Z(I) .LT. 0.0D0) THEN                                       
          Z(I)=Z(I)+RMAX                                                
          Z(I-1)=-1.0D0                                                 
        ENDIF                                                           
278   CONTINUE                                                          
                                                                        
290   IF (Z(1) .GT. 0.5) GOTO 300                                       
      I=1                                                               
291     I=I+1                                                           
        IF (Z(I) .LT. 0.5 .AND. I .LT. L+1) GOTO 291                    
      IF (I .EQ. L+1) THEN                                              
        Z(-1)=1.0D0                                                     
        Z(L+1)=0.0D0                                                    
        GOTO 300                                                        
      ENDIF                                                             
292   DO 293 J=1,L+1-I                                                  
        Z(J)=Z(J+I-1)                                                   
293   CONTINUE                                                          
      DO 294 J=L+2-I,L                                                  
        Z(J)=0.0D0                                                      
294   CONTINUE                                                          
      Z(L+1)=Z(L+1)-I+1                                                 
300   DO 310 I=-1,L+1                                                   
        C(I)=Z(I)                                                       
310   CONTINUE                                                          
311   IF (C(1) .LT. 0.5) THEN                                           
        C(-1)=1.0D0                                                     
        C(L+1)=0.0D0                                                    
      ENDIF                                                             
      RETURN                                                            
END SUBROUTINE ARADD
                                                                        
                                              
SUBROUTINE ARSUB(A,B,C,L,RMAX)                                    
!     ****************************************************************  
!     *                 SUBROUTINE ARSUB                             *  
!     *  Description : Accepts two arrays and subtracts each element *  
!     *    in the second array from the element in the first array   *  
!     *    and returns the solution.  The parameters L and RMAX are  *  
!     *    the size of the array and the number of digits needed for *  
!     *    the accuracy, respectively.                               *  
!     *  Subprograms called: ARADD                                   *  
!     ****************************************************************  
      INTEGER      :: L,I
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)
      REAL(KIND=PR) :: A,B,C,B2,RMAX                                              
      DIMENSION A(-1:*),B(-1:*),C(-1:*),B2(-1:25000)                
                                                                        
      DO 100 I=-1,L+1                                                   
        B2(I)=B(I)                                                      
100   CONTINUE                                                          
      B2(-1)=(-1.0D0)*B2(-1)                                            
      CALL ARADD(A,B2,C,L,RMAX)                                         
      RETURN                                                            
END SUBROUTINE ARSUB
                                                                        
SUBROUTINE ARMULT(A,B,C,L,RMAX)                                   
!     ****************************************************************  
!     *                 SUBROUTINE ARMULT                            *  
!     *  Description : Accepts two arrays and returns the product.   *  
!     *    L and RMAX are the size of the arrays and the number of   *  
!     *    digits needed to represent the numbers with the required  *  
!     *    accuracy.                                                 *  
!     *  Subprograms called: none                                    *  
!     ****************************************************************  
      INTEGER      :: L
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)
      REAL(KIND=PR) :: A,B,C,Z,B2,CARRY,RMAX,RMAX2                                
      DIMENSION A(-1:*),C(-1:*),Z(-1:25000)                           
      INTEGER      :: I                                                         
                                                                        
      RMAX2=1.0D0/RMAX                                                  
      Z(-1)=DSIGN(1.0D0,B)*A(-1)                                        
      B2=DABS(B)                                                        
      Z(L+1)=A(L+1)                                                     
      DO 100 I=0,L                                                      
        Z(I)=0.0D0                                                      
100   CONTINUE                                                          
      IF (B2 .LE. 1.0D-10 .OR. A(1) .LE. 1.0D-10) THEN                  
        Z(-1)=1.0D0                                                     
        Z(L+1)=0.0D0                                                    
        GOTO 198                                                        
      ENDIF                                                             
      DO 110 I=L,1,-1                                                   
        Z(I)=A(I)*B2+Z(I)                                               
        IF (Z(I) .GE. RMAX) THEN                                        
          CARRY=DINT(Z(I)/RMAX)                                         
          Z(I)=Z(I)-CARRY*RMAX                                          
          Z(I-1)=CARRY                                                  
        ENDIF                                                           
110   CONTINUE                                                          
      IF (Z(0) .LT. 0.5) GOTO 150                                       
      DO 120 I=L,1,-1                                                   
        Z(I)=Z(I-1)                                                     
120   CONTINUE                                                          
      Z(L+1)=Z(L+1)+1.0D0                                               
      Z(0)=0.0D0                                                        
150   CONTINUE                                                          
                                                                        
                                                                        
198   DO 199 I=-1,L+1                                                   
        C(I)=Z(I)                                                       
199   CONTINUE                                                          
      IF (C(1) .LT. 0.5) THEN                                           
        C(-1)=1.0D0                                                     
        C(L+1)=0.0D0                                                    
      ENDIF                                                             
      RETURN                                                            
END SUBROUTINE ARMULT
                                                                        
                                                                        
SUBROUTINE CMPADD(AR,AI,BR,BI,CR,CI,L,RMAX)                       
!     ****************************************************************  
!     *                 SUBROUTINE CMPADD                            *  
!     *  Description : Takes two arrays representing one real and    *  
!     *    one imaginary part, and adds two arrays representing      *  
!     *    another complex number and returns two array holding the  *  
!     *    complex sum.                                              *  
!     *              (CR,CI) = (AR+BR, AI+BI)                        *  
!     *  Subprograms called: ARADD                                   *  
!     ****************************************************************  
      INTEGER      :: L
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)

      REAL(KIND=PR) :: AR,AI,BR,BI,CR,CI,RMAX                                     
      DIMENSION AR(-1:*),AI(-1:*),BR(-1:*),BI(-1:*)             
      DIMENSION CR(-1:*),CI(-1:*)                                   
                                                                        
      CALL ARADD(AR,BR,CR,L,RMAX)                                       
      CALL ARADD(AI,BI,CI,L,RMAX)                                       
      RETURN                                                            
END SUBROUTINE CMPADD


SUBROUTINE CMPSUB(AR,AI,BR,BI,CR,CI,L,RMAX)                       
!     ****************************************************************  
!     *                 SUBROUTINE CMPSUB                            *  
!     *  Description : Takes two arrays representing one real and    *  
!     *    one imaginary part, and subtracts two arrays representing *  
!     *    another complex number and returns two array holding the  *  
!     *    complex sum.                                              *  
!     *              (CR,CI) = (AR+BR, AI+BI)                        *  
!     *  Subprograms called: ARADD                                   *  
!     ****************************************************************  
      INTEGER      :: L
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)

      REAL(KIND=PR) :: AR,AI,BR,BI,CR,CI,RMAX                                     
      DIMENSION AR(-1:*),AI(-1:*),BR(-1:*),BI(-1:*)             
      DIMENSION CR(-1:*),CI(-1:*)                                   
                                                                        
      CALL ARSUB(AR,BR,CR,L,RMAX)                                       
      CALL ARSUB(AI,BI,CI,L,RMAX)                                       
      RETURN                                                            
END SUBROUTINE CMPSUB
                                                                        
                                                                        
                                                                        
SUBROUTINE CMPMUL(AR,AI,BR,BI,CR,CI,L,RMAX)                       
!     ****************************************************************  
!     *                 SUBROUTINE CMPMUL                            *  
!     *  Description : Takes two arrays representing one real and    *  
!     *    one imaginary part, and multiplies it with two arrays     *  
!     *    representing another complex number and returns the       *  
!     *    complex product.                                          *  
!     *  Subprograms called: ARMULT, ARSUB, ARADD                    *  
!     ****************************************************************  
      INTEGER      :: L
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)

      REAL(KIND=PR) :: AR,AI,BR,BI,CR,CI,D1,D2,RMAX                           
      DIMENSION AR(-1:*),AI(-1:*),CR(-1:*),CI(-1:*)             
      DIMENSION D1(-1:25000),D2(-1:25000)                       
                                                                        
      CALL ARMULT(AR,BR,D1,L,RMAX)                                      
      CALL ARMULT(AI,BI,D2,L,RMAX)                                      
      CALL ARSUB(D1,D2,CR,L,RMAX)                                      
      CALL ARMULT(AR,BI,D1,L,RMAX)                                      
      CALL ARMULT(AI,BR,D2,L,RMAX)                                      
      CALL ARADD(D1,D2,CI,L,RMAX)                                       
      RETURN                                                            
END SUBROUTINE CMPMUL
                                                                        
                                                                        
                                                                        
SUBROUTINE ARYDIV(AR,AI,BR,BI,C,L,LNCHF,RMAX,BIT)                 
!     ****************************************************************  
!     *                 SUBROUTINE ARYDIV                            *  
!     *  Description : Returns the double precision complex number   *  
!     *    resulting from the division of four arrays, representing  *  
!     *    two complex numbers.  The number returned will be in one  *  
!     *    two different forms.  Either standard scientific or as    *  
!     *    the log of the number.                                    *  
!     *  Subprograms called: CONV21, CONV12, EADD, ECPDIV, EMULT     *  
!     ****************************************************************  
      INTEGER         :: L,BIT,REXP,IR10,II10,LNCHF
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)

      COMPLEX(KIND=PR) :: C                                                      
      REAL(KIND=PR)    :: PHI,N1,N2,N3,E1,E2,E3,RR10,RI10,X              
      REAL(KIND=PR)    :: X1,X2,DUM1,DUM2,RMAX                              
      REAL(KIND=PR), DIMENSION(-1:*) :: AR,AI,BR,BI             
      REAL(KIND=PR), DIMENSION(2,2)  :: AE, BE, CE             
                                                                        
      REXP=BIT/2                                                        
      X=REXP*(AR(L+1)-2)                                                
      RR10=X*DLOG10(2.0D0)/DLOG10(10.0D0)                               
      IR10=INT(RR10)                                                    
      RR10=RR10-IR10                                                    
      X=REXP*(AI(L+1)-2)                                                
      RI10=X*DLOG10(2.0D0)/DLOG10(10.0D0)                               
      II10=INT(RI10)                                                    
      RI10=RI10-II10                                                    
      DUM1=DSIGN(AR(1)*RMAX*RMAX+AR(2)*RMAX+AR(3),AR(-1))               
      DUM2=DSIGN(AI(1)*RMAX*RMAX+AI(2)*RMAX+AI(3),AI(-1))               
      DUM1=DUM1*10**RR10                                                
      DUM2=DUM2*10**RI10                                                
      CALL CONV12(CMPLX(DUM1,DUM2, KIND=PR),AE)                                 
      AE(1,2)=AE(1,2)+IR10                                              
      AE(2,2)=AE(2,2)+II10                                              
      X=REXP*(BR(L+1)-2)                                                
      RR10=X*DLOG10(2.0D0)/DLOG10(10.0D0)                               
      IR10=INT(RR10)                                                    
      RR10=RR10-IR10                                                    
      X=REXP*(BI(L+1)-2)                                                
      RI10=X*DLOG10(2.0D0)/DLOG10(10.0D0)                               
      II10=INT(RI10)                                                    
      RI10=RI10-II10                                                    
      DUM1=DSIGN(BR(1)*RMAX*RMAX+BR(2)*RMAX+BR(3),BR(-1))               
      DUM2=DSIGN(BI(1)*RMAX*RMAX+BI(2)*RMAX+BI(3),BI(-1))               
      DUM1=DUM1*10**RR10                                                
      DUM2=DUM2*10**RI10                                                
      CALL CONV12(CMPLX(DUM1,DUM2, KIND=PR),BE)                                 
      BE(1,2)=BE(1,2)+IR10                                              
      BE(2,2)=BE(2,2)+II10                                              
      CALL ECPDIV(AE,BE,CE)                                             
      IF (LNCHF .EQ. 0) THEN                                            
        CALL CONV21(CE,C)                                               
      ELSE                                                              
        CALL EMULT(CE(1,1),CE(1,2),CE(1,1),CE(1,2),N1,E1)               
        CALL EMULT(CE(2,1),CE(2,2),CE(2,1),CE(2,2),N2,E2)               
        CALL EADD(N1,E1,N2,E2,N3,E3)                                    
        N1=CE(1,1)                                                      
        E1=CE(1,2)-CE(2,2)                                              
        X2=CE(2,1)                                                      
        IF (E1 .GT. 74.0D0) THEN                                        
          X1=1.0D75                                                     
        ELSEIF (E1 .LT. -74.0D0) THEN                                   
          X1=0                                                          
        ELSE                                                            
          X1=N1*(10**E1)                                                
        ENDIF                                                           
        PHI=DATAN2(X2,X1)                                               
        C=CMPLX(0.50D0*(DLOG(N3)+E3*DLOG(10.0D0)),PHI)                 
      ENDIF                                                             
      RETURN                                                            
END SUBROUTINE ARYDIV


SUBROUTINE EMULT(N1,E1,N2,E2,NF,EF)                               
!     ****************************************************************  
!     *                 SUBROUTINE EMULT                             *  
!     *  Description : Takes one base and exponent and multiplies it *  
!     *    by another numbers base and exponent to give the product  *  
!     *    in the form of base and exponent.                         *  
!     *  Subprograms called: none                                    *  
!     ****************************************************************  
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)
      REAL(KIND=PR) :: N1,E1,N2,E2,NF,EF                                          
                                                                        
      NF=N1*N2                                                          
      EF=E1+E2                                                          
      IF (DABS(NF) .GE. 10.0D0) THEN                                    
        NF=NF/10.0D0                                                    
        EF=EF+1.0D0                                                     
      ENDIF                                                             
      RETURN                                                            
END SUBROUTINE EMULT


                                                                        
SUBROUTINE EDIV(N1,E1,N2,E2,NF,EF)                                
!     ****************************************************************  
!     *                 SUBROUTINE EDIV                              *  
!     *  Description : returns the solution in the form of base and  *  
!     *    exponent of the division of two exponential numbers.      *  
!     *  Subprograms called: none                                    *  
!     ****************************************************************  
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)
      REAL(KIND=PR) :: N1,E1,N2,E2,NF,EF                                          
                                                                        
      NF=N1/N2                                                          
      EF=E1-E2                                                          
      IF ((DABS(NF) .LT. 1.0D0) .AND. (NF .NE. 0.0D0)) THEN             
        NF=NF*10.0D0                                                    
        EF=EF-1.0D0                                                     
      ENDIF                                                             
      RETURN                                                            
END SUBROUTINE EDIV


SUBROUTINE EADD(N1,E1,N2,E2,NF,EF)                                
!     ****************************************************************  
!     *                 SUBROUTINE EADD                              *  
!     *  Description : Returns the sum of two numbers in the form    *  
!     *    of a base and an exponent.                                *  
!     *  Subprograms called: none                                    *  
!     ****************************************************************  
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)
      REAL(KIND=PR) :: N1,E1,N2,E2,NF,EF,EDIFF                                    

      EDIFF=E1-E2                                                       
      IF (EDIFF .GT. 36.0D0) THEN                                       
        NF=N1                                                           
        EF=E1                                                           
      ELSE IF (EDIFF .LT. -36.0D0) THEN                                 
        NF=N2                                                           
        EF=E2                                                           
      ELSE                                                              
        NF=N1*(10.0D0**EDIFF)+N2                                        
        EF=E2                                                           
400     IF (DABS(NF) .LT. 10.0D0) GOTO 410                              
          NF=NF/10.0D0                                                  
          EF=EF+1.0D0                                                   
          GOTO 400                                                      
410     IF ((DABS(NF) .GE. 1.0D0) .OR. (NF .EQ. 0.0D0)) GOTO 420        
          NF=NF*10.0D0                                                  
          EF=EF-1.0D0                                                   
          GOTO 410                                                      
      ENDIF                                                             
420   RETURN                                                            
END SUBROUTINE EADD


SUBROUTINE ESUB(N1,E1,N2,E2,NF,EF)                                
!     ****************************************************************  
!     *                 SUBROUTINE ESUB                              *  
!     *  Description : Returns the solution to the subtraction of    *  
!     *    two numbers in the form of base and exponent.             *  
!     *  Subprograms called: EADD                                    *  
!     ****************************************************************  
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)
      REAL(KIND=PR) :: N1,E1,N2,E2,NF,EF                                          
                                                                        
      CALL EADD(N1,E1,N2*(-1.0D0),E2,NF,EF)                             
      RETURN                                                            
END SUBROUTINE ESUB


SUBROUTINE CONV12(CN,CAE)                                         
!     ****************************************************************  
!     *                 SUBROUTINE CONV12                            *  
!     *  Description : Converts a number from complex notation to a  *  
!     *    form of a 2x2 real array.                                 *  
!     *  Subprograms called: none                                    *  
!     ****************************************************************  
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)
      COMPLEX(KIND=PR) :: CN                                                     
      REAL(KIND=PR), DIMENSION(2,2)    :: CAE
                                                                        
      CAE(1,1)=DBLE(CN)                                                
      CAE(1,2)=0.0D0                                                    
300   IF (DABS(CAE(1,1)) .LT. 10.0D0) GOTO 310                          
        CAE(1,1)=CAE(1,1)/10.0D0                                        
        CAE(1,2)=CAE(1,2)+1.0D0                                         
        GOTO 300                                                        
310   IF ((DABS(CAE(1,1)) .GE. 1.0D0) .OR. (CAE(1,1) .EQ. 0.0D0)) GOTO 320                                                         
        CAE(1,1)=CAE(1,1)*10.0D0                                        
        CAE(1,2)=CAE(1,2)-1.0D0                                         
        GOTO 310                                                        
320   CAE(2,1)=AIMAG(CN)                                                
      CAE(2,2)=0.0D0                                                    
330   IF (DABS(CAE(2,1)) .LT. 10.0D0) GOTO 340                          
        CAE(2,1)=CAE(2,1)/10.0D0                                        
        CAE(2,2)=CAE(2,2)+1.0D0                                         
        GOTO 330                                                        
340   IF ((DABS(CAE(2,1)) .GE. 1.0D0) .OR. (CAE(2,1) .EQ. 0.0D0)) GOTO 350                                                         
        CAE(2,1)=CAE(2,1)*10.0D0                                        
        CAE(2,2)=CAE(2,2)-1.0D0                                         
        GOTO 340                                                        
350   RETURN                                                            
END SUBROUTINE CONV12


SUBROUTINE CONV21(CAE,CN)                                         
!     ****************************************************************  
!     *                 SUBROUTINE CONV21                            *  
!     *  Description : Converts a number represented in a 2x2 real   *  
!     *    array to the form of a complex number.                    *  
!     *  Subprograms called: none                                    *  
!     ****************************************************************  
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)
      REAL(KIND=PR), DIMENSION(2,2)     :: CAE           
      COMPLEX(KIND=PR) :: CN                                                     

      IF (CAE(1,2) .GT. 75 .OR. CAE(2,2) .GT. 75) THEN                  
        CN=CMPLX(1.0D75,1.0D75, KIND=PR)                                         
      ELSE IF (CAE(2,2) .LT. -75) THEN                                  
        CN=CMPLX(CAE(1,1)*(10**CAE(1,2)),0D0)                          
      ELSE                                                              
        CN=CMPLX(CAE(1,1)*(10**CAE(1,2)),CAE(2,1)*(10**CAE(2,2)))      
        CN=CMPLX(CAE(1,1)*(10**CAE(1,2)),CAE(2,1)*(10**CAE(2,2)))      
      ENDIF                                                             
      RETURN                                                            
END SUBROUTINE CONV21


SUBROUTINE ECPMUL(A,B,C)                                          
!     ****************************************************************  
!     *                 SUBROUTINE ECPMUL                            *  
!     *  Description : Multiplies two numbers which are each         *  
!     *    represented in the form of a two by two array and returns *  
!     *    the solution in the same form.                            *  
!     *  Subprograms called: EMULT, ESUB, EADD                       *  
!     ****************************************************************  
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)
      REAL(KIND=PR)  :: N1,E1,N2,E2                                    
      REAL(KIND=PR), DIMENSION(2,2) :: A,B,C,C2                            
                                                                        
      CALL EMULT(A(1,1),A(1,2),B(1,1),B(1,2),N1,E1)                     
      CALL EMULT(A(2,1),A(2,2),B(2,1),B(2,2),N2,E2)                     
      CALL ESUB(N1,E1,N2,E2,C2(1,1),C2(1,2))                            
      CALL EMULT(A(1,1),A(1,2),B(2,1),B(2,2),N1,E1)                     
      CALL EMULT(A(2,1),A(2,2),B(1,1),B(1,2),N2,E2)                     
      CALL EADD(N1,E1,N2,E2,C(2,1),C(2,2))                              
      C(1,1)=C2(1,1)                                                    
      C(1,2)=C2(1,2)                                                    
      RETURN                                                            
END SUBROUTINE ECPMUL


SUBROUTINE ECPDIV(A,B,C)                                          
!     ****************************************************************  
!     *                 SUBROUTINE ECPDIV                            *  
!     *  Description : Divides two numbers and returns the solution. *  
!     *    All numbers are represented by a 2x2 array.               *  
!     *  Subprograms called: EADD, ECPMUL, EDIV, EMULT               *  
!     ****************************************************************  
      INTEGER, PARAMETER    :: PR=KIND(1.0D0)
      REAL(KIND=PR)  ::  N1,E1,N2,E2,N3,E3
      REAL(KIND=PR), DIMENSION(2,2)  ::  A,B,C,B2,C2                              
                                                                        
      B2(1,1)=B(1,1)                                                    
      B2(1,2)=B(1,2)                                                    
      B2(2,1)=-1.0D0*B(2,1)                                             
      B2(2,2)=B(2,2)                                                    
      CALL ECPMUL(A,B2,C2)                                              
      CALL EMULT(B(1,1),B(1,2),B(1,1),B(1,2),N1,E1)                     
      CALL EMULT(B(2,1),B(2,2),B(2,1),B(2,2),N2,E2)                     
      CALL EADD(N1,E1,N2,E2,N3,E3)                                      
      CALL EDIV(C2(1,1),C2(1,2),N3,E3,C(1,1),C(1,2))                    
      CALL EDIV(C2(2,1),C2(2,2),N3,E3,C(2,1),C(2,2))                    
      RETURN                                                            
END SUBROUTINE ECPDIV

                                                             

