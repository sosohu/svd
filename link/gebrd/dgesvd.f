      PROGRAM TEST_DGESVD
*     .. Scalar Arguments ..
*     ..
*  =====================================================================
*
*
*
*     .. Parameters ..  
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER             M,N,LDA
      PARAMETER          (M=400,N=400,LDA=M+1)
*     ..
*     .. Local Scalars ..
      REAL               START,FINISH,XX
      INTEGER            I,J,K,L,TIMES,INFO,LWORK,LDU,LDVT
      PARAMETER          ( MINMN = MIN(M,N),LDU=M,LDVT=N,
     $              LWORK = MAX(3*MIN(M,N)+MAX(M,N),5*MIN(M,N)) )
      DOUBLE PRECISION   RUNTIMES
      DOUBLE PRECISION   UVT
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   A(LDA,N)
      DOUBLE PRECISION   S(MINMN),U(LDU,M),VT(LDVT,N),WORK(LWORK)
*     ..
*     .. External Functions ..
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGESVD
*     ..
*     .. Scalars in Common ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
      TIMES=1
      call random_seed()
      DO I = 1,M
         DO J = 1,N
            call random_number(XX)
            A(I,J)=XX
         END DO
      END DO
         
      call cpu_time(start)
      CALL DGESVD( 'A','A',M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,INFO)
      call cpu_time(finish)
      RUNTIMES=(finish-start)*1000000
      WRITE(*,200)N,TIMES,RUNTIMES,INFO

	  DO I = 1,MINMN
           WRITE(*,"(F10.6)",advance='no') S(I)
	  END DO
 
	  PRINT*, ' '
	  PRINT*, 'U*UT'

	  DO I = 1,M
	   DO J = 1,M
	       UVT = 0
	       DO K = 1,M
           UVT = UVT + U(I,K)*U(J,K)  
		   END DO
           WRITE(*,"(F12.9)",advance='no') UVT
	   END DO
       WRITE(*,*) ' '
	  END DO

	  PRINT*, 'V*VT'

	  DO I = 1,N
	   DO J = 1,N
	       UVT = 0
		   DO K = 1,N
		   UVT = UVT + VT(I,K)*VT(J,K)
		   END DO
           WRITE(*,"(F12.9)",advance='no') UVT
	   END DO
       WRITE(*,*) ' '
	  END DO

      RETURN

 200  FORMAT (' N=',I4,' TIMES=',I5,' RUNTIMES=',F20.2,' INFO=',I1)
*
*     End of TEST_DGESVD
*
      END
