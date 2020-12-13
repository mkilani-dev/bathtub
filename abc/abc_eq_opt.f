      PROGRAM ABC
      IMPLICIT NONE

      INTEGER LDA
      PARAMETER (LDA=1000)
C Model parameters and functions
      DOUBLE PRECISION L, VF, CAP, AA, BB, GG
      DOUBLE PRECISION T, DT
C Values for DNSQE routine
      DOUBLE PRECISION SMCOST(LDA), K(LDA), JAK(LDA,LDA),Z(LDA)
      INTEGER N, C, IT, IPVT(LDA), JOB
      DOUBLE PRECISION ZERO, ONE, X, SMC, JSMC, B(LDA),RCOND,SC,KK
      DOUBLE PRECISION SOMT, SOMK, CST, UC

      INTEGER I, J
      COMMON /PARAM/ L, VF, CAP, AA, BB, GG
      COMMON /CENTRAL/ C

      DATA L,VF,CAP,AA,BB,GG /5.D0,15.D0,6.D0,10.D0,8.D0,15.D0/


      PRINT*, "============================="
      PRINT*, "====== OPTIMUM =============="
      PRINT*, "============================="

      CALL ONEMASS(X,SC)
      PRINT '(2F9.3)', X, SC 

      N = 2
      C = 2

      K(1) = 0.0D0
      K(2) = X

      CALL DISPLAY(N,C,K)
      DO 10, I = 1, 3
      CALL ITERATION(N,C,K,SC) 
      CALL DISPLAY(N,C,K)
 10   CONTINUE 

C     KK = 2.205D0
C     PRINT*, N
C     KK = K(4)
C     PRINT*, AA*(KK*DT(KK)+T(KK)) + GG*SOMT(N,K,4,N-1)+
C    &                               GG*DT(KK)*SOMK(N,K,4,N-1)
C     PRINT*, C
C     PRINT*, SMC(N,C,K,4)
C     PRINT*, K(1:4)
      PRINT*, "============================="
      PRINT*, "====== EQUILIBRIUM =========="
      PRINT*, "============================="
      CALL EQONEMASS(K,CST,C)
      CALL EQDISPLAY(2,C,K)

      N = 2
C     CALL EQCSOLVE(N,C,K,8.5D0)
      CALL EQITERATION(N,C,K,SC) 
      CALL EQDISPLAY(N,C,K)
      PRINT*, "== OK =="

      CALL EQITER2(N,C,K,SC)
C     PRINT*, SC
C     CALL EQCSOLVE(N,C,K,49.13D0)
C     CALL EQDISPLAY(N,C,K)
C     PRINT*, AA*T(0.D0)+BB*(T(K(1))+T(K(2)) )
      CALL EQITERATION(N,C,K,SC) 
      CALL EQDISPLAY(N,C,K)
      CALL EQCSOLVE(N, C, K, 35.0D0) 
      CALL EQDISPLAY(N,C,K)
      CALL EQCSOLVE(N, C, K, 45.0D0) 
      CALL EQDISPLAY(N,C,K)
      CALL EQCSOLVE(N, C, K, 52.075D0) 
      CALL EQDISPLAY(N,C,K)
      CALL EQCSOLVE(N, C, K, 52.0833D0) 
      CALL EQDISPLAY(N,C,K)

      N = N + 1
      K(N) = 0.D0
      CALL EQCSOLVE(N, C, K, 83.3033D0) 
      CALL EQDISPLAY(N,C,K)


      N = N + 1
      K(7) = 0.D0
      K(6) = K(5)
      K(5) = K(4)
      K(4) = K(3)
      K(3) = K(2)
      K(2) = K(1)
      K(1) = 0.D0
      C = C + 1
      CALL EQDISPLAY(N,C,K)
      CALL EQCSOLVE(N, C, K, 100.6833D0) 
      CALL EQDISPLAY(N,C,K)

C     PRINT*, UC(N,C,K,5)-UC(N,C,K,4)

C     CALL EQITERATION(N,C,K,SC) 
C     CALL EQITER2(N,C,K,SC)
C     CALL EQDISPLAY(N,C,K)


C     PRINT*, "UC ---------> ", SC
C     CALL EQITER2(N,C,K,SC) 
C     PRINT*, "UC ---------> ", SC
C     CALL EQDISPLAY(N,C,K)
C     N = N + 1
C     K(N) = 0.0D0
C     CALL EQCSOLVE(N, C, K, 55.075D0) 
C     CALL EQDISPLAY(N,C,K)
C     CALL EQCSOLVE(N, C, K, 85.075D0) 
C     CALL EQDISPLAY(N,C,K)
C     CALL EQCSOLVE(N, C, K, SC) 
C     PRINT*, AA*T(0.0D0)+GG*(T(K(3))+T(K(4)))

C     CALL EQITERATION(N,C,K,SC) 
C     CALL EQDISPLAY(N,C,K)

      END 

      SUBROUTINE DISPLAY (N,C,K)
      IMPLICIT NONE
      INTEGER N, C, I
      DOUBLE PRECISION K(N), SMC, T
      DOUBLE PRECISION SOMK

      PRINT*
      PRINT '(I5,A1,2F9.3)', 0,' ', 0.00D0 , SMC(N,C,K,0)
      DO 10, I = 1, N 
      IF(I.EQ.C) THEN
        PRINT '(I5,A1,3F9.3)', I,'*', K(I), SMC(N,C,K,I), T(K(I)) 
      ELSE
        PRINT '(I5,A1,3F9.3)', I,' ', K(I), SMC(N,C,K,I), T(K(I))
      ENDIF
 10   CONTINUE
      PRINT '(I5,A1,2F9.3)', N+1,' ', 0.00D0 , SMC(N,C,K,N+1)
      PRINT*
      PRINT '(A12,F9.3)', "Population :", SOMK(N,K,1,N)

      RETURN
      END

      SUBROUTINE EQDISPLAY (N,C,K)
      IMPLICIT NONE
      INTEGER N, C, I
      DOUBLE PRECISION K(N), UC, T
      DOUBLE PRECISION SOMK

      PRINT*
      PRINT '(I5,A1,2F9.3)', 0,' ', 0.00D0 , UC(N,C,K,0)
      DO 10, I = 1, N 
      IF(I.EQ.C) THEN
        PRINT '(I5,A1,3F9.3)', I,'*', K(I), UC(N,C,K,I), T(K(I)) 
      ELSE
        PRINT '(I5,A1,3F9.3)', I,' ', K(I), UC(N,C,K,I), T(K(I))
      ENDIF
 10   CONTINUE
      PRINT '(I5,A1,2F9.3)', N+1,' ', 0.00D0 , UC(N,C,K,N+1)
      PRINT*
      PRINT '(A12,F9.3)', "Population :", SOMK(N,K,1,N)

      RETURN
      END

      SUBROUTINE ITERATION(N,C,K,SC)
      IMPLICIT NONE
      INTEGER N, C, I
      DOUBLE PRECISION SC, VAL, K(N+1), SMC

      CALL ITER1(N,C,K,SC)
      CALL CSOLVE(N, C, K, SC) 
      VAL = SMC(N,C,K,N) - SMC(N,C,K,N+1) 
      IF (VAL.LT.0.0D0) THEN
        DO 3, I = N+1, 2, -1
 3      K(I) = K(I-1)
        K(1) = 0.0D0
        C = C + 1
        N = N + 1
      ELSE
        CALL ITER2(N,C,K,SC)
        CALL CSOLVE(N, C, K, SC) 
        N = N + 1
        K(N) = 0.0D0
      ENDIF 

      RETURN
      END

      SUBROUTINE ITER1(N,C,K,SC)
C Find marginal social cost SC such that mass
C 0 (late arrival) is just to become active
C Newton algorithm is used for this one dimensional problem. 
C The jacobian is approximated
C using finite differences.
      IMPLICIT NONE
      INTEGER N, C, LDA, IT,ITMAX
      PARAMETER(LDA=1000)
      DOUBLE PRECISION SMC, K(N), H, VAL0,VAL1,VAL2,X,DIFF,SC,TOL
      DOUBLE PRECISION TMP

      ITMAX = 25

      H = 0.01D0
      X = SC
      IT = 0
      TOL = 1.0D-9
  
  1   CONTINUE

      IT = IT + 1

      TMP = SMC(N,C,K,0)

      CALL CSOLVE(N,C,K,X)
      VAL0 = X - SMC(N,C,K,0)

      CALL CSOLVE(N,C,K,X+H)
      VAL1 = X+H - SMC(N,C,K,0)
      CALL CSOLVE(N,C,K,X-H)
      VAL2 = X-H - SMC(N,C,K,0)

      DIFF = (VAL2-VAL1)/(2.0D0*H)

      X = X + VAL0 / DIFF 

C     PRINT *, X, VAL0, DIFF, IT, K(1:3)
      IF(ABS(VAL0).GT.TOL.AND.IT.LT.ITMAX) GOTO 1 

      SC = X

      RETURN
      END

      SUBROUTINE ITER2(N,C,K,SC)
C Find marginal social cost SC such that mass
C N+1 (late arrival) is just to become active
C Newton algorithm is used for this one dimensional problem. 
C The jacobian is approximated
C using finite differences.
      IMPLICIT NONE
      INTEGER N, C, LDA, IT,ITMAX
      PARAMETER(LDA=1000)
      DOUBLE PRECISION SMC, K(N), H, VAL0,VAL1,VAL2,X,DIFF,SC,TOL
      DOUBLE PRECISION TMP

      ITMAX = 25

      H = 0.01D0
      X = SC
      IT = 0
      TOL = 1.0D-9
  
  1   CONTINUE

      IT = IT + 1

      TMP = SMC(N,C,K,N+1)

      CALL CSOLVE(N,C,K,X)
      VAL0 = X - SMC(N,C,K,N+1)

      CALL CSOLVE(N,C,K,X+H)
      VAL1 = X+H - SMC(N,C,K,N+1)
      CALL CSOLVE(N,C,K,X-H)
      VAL2 = X-H - SMC(N,C,K,N+1)

      DIFF = (VAL2-VAL1)/(2.0D0*H)

      X = X + VAL0 / DIFF 

C     PRINT *, X, VAL0, DIFF, IT, K(1:3)
      IF(ABS(VAL0).GT.TOL.AND.IT.LT.ITMAX) GOTO 1 

      SC = X

      RETURN
      END

      SUBROUTINE CSOLVE(N,C,K,SC)
C
C Given N and C find K(1)..K(N) such that
C the marginal social cost is equal to SC;
C The computation is performed with active 
C being unchanged.
C 
      IMPLICIT NONE 
      INTEGER N, C, LDA, I, J, JOB, IT
      PARAMETER (LDA=1000)
      DOUBLE PRECISION K(LDA), JAK(LDA,LDA),SMCOST(LDA), TOL
      DOUBLE PRECISION B(LDA), RCOND, Z(N), SC, SMC, JSMC, DENORM
      INTEGER IPVT(N), ITMAX

      ITMAX = 25
      TOL = 1.0D-9

      IT = 0

1     CONTINUE

      IT = IT + 1

      DO 10, I = 1, N
      B(I) = SMC(N,C,K,I) - SC
      DO 10, J = 1, N
      JAK(I,J) = JSMC(N,C,K,I,J)
10    CONTINUE

      JOB = 0
      CALL DGECO(JAK,LDA,N,IPVT,RCOND,Z)
      CALL DGESL(JAK,LDA,N,IPVT,B,JOB) 

      DO 20, I = 1, N
      K(I) = K(I) - B(I)
C     PRINT*, IT, I, K(I), B(I), DENORM(N, B)
 20   CONTINUE

      IF(DENORM(N,B).LT.TOL) GOTO 30

      IF(IT.LT.ITMAX) GOTO 1

      PRINT*, "Maximum number of iterations reached !"

 30   CONTINUE

      RETURN
      END

      SUBROUTINE SMC1(N,C,K,SMCOST)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION K(N), SMCOST(N+2), T, DT
      INTEGER I, J, C
      DOUBLE PRECISION L, VF, CAP, AA, BB, GG
      DOUBLE PRECISION TMP1, TMP2, TMP3, ZERO, SOMT, SOMK

      COMMON /PARAM/   L, VF, CAP, AA, BB, GG

      ZERO = 0.0D0

C MASSES THAT ARRIVE EARLY
      DO 10, I = 1, C-1

      SMCOST(I) = AA * ( K(I) * DT( K(I) ) + T( K(I)) )
     &        + BB *            SOMT(N,K,I+1,C)
     &        + BB * DT(K(I)) * SOMK(N,K,1,I-1) 
     &        - K(N)

 10   CONTINUE

C THE MASS THAT ARRIVES ON TIME (MASS C)

      SMCOST(C) = AA * ( K(C) * DT( K(C) ) + T( K( C )) )  
     &        + BB * DT(K(C)) * SOMK(N,K,1,C-1)
     &        - K(N)


C THE MASSES THAT ARRIVE LATE
      DO 30, I = C + 1, N  

      SMCOST(I) = AA * ( K(I) * DT( K(I) ) + T( K(I)) ) 
     &        + GG * SOMT(N,K,C+1,I)
     &        + GG * DT(K(I)) * SOMK(N,K,I,N-1)
     &        - K(N)

 30   CONTINUE

C THE NEW EARLY MASS 

      SMCOST(N+1) = AA * T(ZERO) + BB * SOMT(N,K,1,C) - K(N)

C THE NEW LATE MASS 

      SMCOST(N+2) = AA * T(ZERO) + GG * (SOMT(N,K,C+1,N)+T(ZERO)) - K(N)


      RETURN
      END

      FUNCTION F(K)
      IMPLICIT NONE

      DOUBLE PRECISION K, F
      DOUBLE PRECISION T, DT
      DOUBLE PRECISION L, VF, CAP, AA, BB, GG

      COMMON /PARAM/ L, VF, CAP, AA, BB, GG 

      F = AA * ( K * DT(K) + T(K) - T(0.0D0) ) - BB * T(K)

      RETURN
      END

      FUNCTION T(K)
      IMPLICIT NONE

      DOUBLE PRECISION K
      DOUBLE PRECISION T
      DOUBLE PRECISION L, VF, CAP, AA, BB, GG

      COMMON /PARAM/ L, VF, CAP, AA, BB, GG

      T = L / ( VF * ( 1.0D0 - K / CAP ) )

      RETURN
      END

      FUNCTION DT(K)
      IMPLICIT NONE

      DOUBLE PRECISION K
      DOUBLE PRECISION DT, T
      DOUBLE PRECISION L, VF, CAP, AA, BB, GG

      COMMON /PARAM/ L, VF, CAP, AA, BB, GG

      DT = VF / ( CAP * L ) * T ( K ) ** 2

      RETURN
      END

      FUNCTION DDT(K)
      IMPLICIT NONE

      DOUBLE PRECISION K
      DOUBLE PRECISION DDT, DT, T
      DOUBLE PRECISION L, VF, CAP, AA, BB, GG

      COMMON /PARAM/ L, VF, CAP, AA, BB, GG

      DDT = 2.0D0 * VF / ( CAP * L ) * T ( K ) * DT ( K )

      RETURN
      END

      SUBROUTINE ONEMASS(X, CST)
      IMPLICIT NONE
      DOUBLE PRECISION F, X, CST
      DOUBLE PRECISION B, C, R, RE, AE
      DOUBLE PRECISION L, VF, CAP, AA, BB, GG
      DOUBLE PRECISION T, DT
      INTEGER IFLAG
      EXTERNAL F
      COMMON /PARAM/ L, VF, CAP, AA, BB, GG

      B = 0.D0
      C = CAP * 0.9D0
      R = CAP / 2.0D0
      RE = 1.0D-9
      AE = 1.0D-9

      CALL DFZERO(F, B, C, R, RE, AE, IFLAG)
      CST = AA*( B * DT( B ) + T( B ) )

      X = B

      RETURN
      END

      FUNCTION SOMT(N,K,I0,I1)
      IMPLICIT NONE
      INTEGER N, I0, I1, I
      DOUBLE PRECISION K(N), TMP, SOMT, T

      TMP = 0.0D0
      DO 1, I = I0, I1
      TMP = TMP + T(K(I))
 1    CONTINUE

      SOMT = TMP

      RETURN
      END

      FUNCTION SOMK(N,K,I0,I1)
      IMPLICIT NONE
      INTEGER N, I0, I1, I
      DOUBLE PRECISION K(N), TMP, SOMK

      TMP = 0.0D0
      DO 1, I = I0, I1
      TMP = TMP + K(I)
 1    CONTINUE

      SOMK = TMP

      RETURN
      END

      FUNCTION JSMC(N,C,K,I,J)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION K(N), K1(N),K2(N), T, DT
      INTEGER I, J, C, IT
      DOUBLE PRECISION L, VF, CAP, AA, BB, GG
      DOUBLE PRECISION SOMT, SOMK, SMC, JSMC, H

      COMMON /PARAM/   L, VF, CAP, AA, BB, GG

      H = 0.01D0

      DO 1, IT = 1, N
      K1(IT) = K(IT)
      K2(IT) = K(IT)
  1   CONTINUE
      K1(J) = K(J) + H
      K2(J) = K(J) - H

      JSMC = ( SMC(N,C,K1,I) - SMC(N,C,K2,I) ) / (2.0D0 * H)

      RETURN
      END

      FUNCTION SMC(N,C,K,I)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION K(N), T, DT
      INTEGER I, C
      DOUBLE PRECISION L, VF, CAP, AA, BB, GG
      DOUBLE PRECISION SOMT, SOMK, SMC

      COMMON /PARAM/   L, VF, CAP, AA, BB, GG

      IF (I.GE.1.AND.I.LT.C) THEN
         SMC = AA * ( K(I) * DT( K(I) ) + T( K(I)) )
     &       + BB *            SOMT(N,K,I+1,C)
     &       + BB * DT(K(I)) * SOMK(N,K,1,I-1) 
      ELSEIF(I.EQ.C) THEN 
         SMC  = AA * ( K(C) * DT( K(C) ) + T( K( C )) )  
     &        + BB * DT(K(C)) * SOMK(N,K,1,C-1)
      ELSEIF(I.GT.C.AND.I.LE.N) THEN
         SMC  = AA * ( K(I) * DT( K(I) ) + T( K(I)) ) 
     &        + GG *            SOMT(N,K,C+1,I)
     &        + GG * DT(K(I)) * SOMK(N,K,I,N)
      ELSEIF(I.EQ.0) THEN
         SMC  = AA * T(0.0D0) + BB * SOMT(N,K,1,C)
      ELSEIF(I.EQ.N+1) THEN
         SMC  = AA * T(0.0D0) + GG * SOMT(N,K,C+1,N+1)
      ELSE
         SMC = 0.0D0
         PRINT*, "Inconsistent call to SCM"
      ENDIF

      RETURN
      END

      FUNCTION UC(N,C,K,I)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION K(N), T, DT
      INTEGER I, C
      DOUBLE PRECISION L, VF, CAP, AA, BB, GG
      DOUBLE PRECISION SOMT, SOMK, UC

      COMMON /PARAM/   L, VF, CAP, AA, BB, GG

      IF (I.GE.1.AND.I.LT.C) THEN
         UC  = AA * T( K(I)) + BB * SOMT(N,K,I+1,C)
      ELSEIF(I.EQ.C) THEN 
         UC  = AA * T( K( C ) )   
      ELSEIF(I.GT.C.AND.I.LE.N) THEN
         UC  = AA * T( K(I) ) + GG * SOMT(N,K,C+1,I)
      ELSEIF(I.EQ.0) THEN
         UC  = AA * T(0.0D0) + BB * SOMT(N,K,1,C)
      ELSEIF(I.EQ.N+1) THEN
         UC  = AA * T(0.0D0) + GG * SOMT(N,K,C+1,N+1)
      ELSE
         UC  = 0.0D0
         PRINT*, "Inconsistent call to UC"
      ENDIF

      RETURN
      END

      FUNCTION JUC(N,C,K,I,J)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION K(N), K1(N),K2(N), T, DT
      INTEGER I, J, C, IT
      DOUBLE PRECISION L, VF, CAP, AA, BB, GG
      DOUBLE PRECISION SOMT, SOMK, SMC, JUC, UC, H

      COMMON /PARAM/   L, VF, CAP, AA, BB, GG

      H = 0.01D0

      DO 1, IT = 1, N
      K1(IT) = K(IT)
      K2(IT) = K(IT)
  1   CONTINUE
      K1(J) = K(J) + H
      K2(J) = K(J) - H

      JUC = ( UC(N,C,K1,I) - UC(N,C,K2,I) ) / (2.0D0 * H)

      RETURN
      END

      SUBROUTINE EQONEMASS(K, CST, C)
C Compute the equilibrium value of oppulatio where the second
C mass after 0 (the central mass) becomes active
      IMPLICIT NONE
      DOUBLE PRECISION F, X, CST, K(*)
      DOUBLE PRECISION L, VF, CAP, AA, BB, GG
      DOUBLE PRECISION T, DT, DIFF, TOL
      INTEGER IFLAG, IT, ITMAX, C
      EXTERNAL F
      COMMON /PARAM/ L, VF, CAP, AA, BB, GG

      TOL=1.0D-9
      ITMAX = 25

C initial solution 
      X = CAP/2.0D0 
C
C 1. Use Newton iterates to find the zero of equation :
C    ALPHA * T(K_0) =  ALPHA T(0) + BETA T(K_0)
C or
C    T(K_0) - ALPHA/(ALPHA-BETA) T(0) = 0
C
C 2. Use Newton iterates to find the zero of equation :
C    ALPHA  T(K_0) =  (ALPHA + GAMMA)  T(0) 
C or
C    T(K_0) - (ALPHA+GAMMA)/ALPHA T(0) = 0
C
      IT = 0
      IF(GG.GE.AA*BB/(AA-BB)) GOTO 10
      IF(GG.LT.AA*BB/(AA-BB)) GOTO 20

  10  CONTINUE
      IT = IT + 1
      DIFF = ( T(X) - AA * T(0.0D0) / (AA-BB) ) / DT(X)
      X = X - DIFF
      IF(ABS(DIFF).LT.TOL  ) GOTO 11
      IF(IT  .LT.ITMAX) GOTO 10
      PRINT*, "Maximum number of iterations reached by EQONEMASS"
  11  CONTINUE
      K(1) = 0.0D0
      K(2) = X
      C = 2 
      GOTO 100 

  20  CONTINUE
      IT = IT + 1
      DIFF = ( T(X) - (AA+GG) * T(0.0D0) / AA ) / DT(X)
      X = X - DIFF
      IF(ABS(DIFF).LT.TOL  ) GOTO 21
      IF(     IT  .LT.ITMAX) GOTO 20
      PRINT*, "Maximum number of iterations reached by EQONEMASS"
  21  CONTINUE
      K(1) = X
      K(2) = 0.0D0
      C = 1
      GOTO 100 

 100  CONTINUE

      CST = AA * T(X)

      RETURN
      END 

      SUBROUTINE EQITERATION(N,C,K,SC)
      IMPLICIT NONE
      INTEGER N, C, I
      DOUBLE PRECISION SC, VAL, K(N+1), UC

      CALL EQITER2(N,C,K,SC)
      CALL EQCSOLVE(N, C, K, SC) 
      VAL = UC(N,C,K,1) - UC(N,C,K,0) 
      IF (VAL.LT.0.0D0) THEN
        N = N + 1
        K(N+1) = 0.0D0
      ELSE
        CALL EQITER1(N,C,K,SC)
        CALL EQCSOLVE(N, C, K, SC) 
        DO 3, I = N+1, 2, -1
 3      K(I) = K(I-1)
        K(1) = 0.0D0
        C = C + 1
        N = N + 1
      ENDIF 

      RETURN
      END

      SUBROUTINE EQCSOLVE(N,C,K,CBAR)
C
C Given N and C find K(1)..K(N) such that
C the marginal social cost is equal to CBAR;
C The computation is performed with the 
C set of active masses being unchanged.
C 
      IMPLICIT NONE 
      INTEGER N, C, LDA, I, J, JOB, IT
      PARAMETER (LDA=1000)
      DOUBLE PRECISION K(LDA), JAK(LDA,LDA), TOL
      DOUBLE PRECISION B(LDA), RCOND, Z(N), CBAR, UC, JUC, DENORM
      INTEGER IPVT(N), ITMAX

      ITMAX = 25
      TOL = 1.0D-9

      IT = 0

1     CONTINUE

      IT = IT + 1

      DO 10, I = 1, N
      B(I) = UC(N,C,K,I) - CBAR
      DO 10, J = 1, N
      JAK(I,J) = JUC(N,C,K,I,J)
10    CONTINUE

      JOB = 0
      CALL DGECO(JAK,LDA,N,IPVT,RCOND,Z)
      CALL DGESL(JAK,LDA,N,IPVT,B,JOB) 

      DO 20, I = 1, N
      K(I) = K(I) - B(I)
C     PRINT*, IT, I, K(I), B(I), DENORM(N, B), CBAR
 20   CONTINUE

      IF(DENORM(N,B).LT.TOL) GOTO 30

      IF(IT.LT.ITMAX) GOTO 1

      PRINT*, "Maximum number of iterations reached !"

 30   CONTINUE

      RETURN
      END

      SUBROUTINE EQITER1(N,C,K,SC)
C Find marginal social cost SC such that mass
C 0 (early arrival) is just to become active
C Newton algorithm is used for this one dimensional 
C problem. 
C The jacobian is approximated
C using finite differences.
      IMPLICIT NONE
      INTEGER N, C, LDA, IT,ITMAX
      PARAMETER(LDA=1000)
      DOUBLE PRECISION UC, K(N), H, VAL0,VAL1,VAL2,X,DIFF,SC,TOL
      DOUBLE PRECISION TMP

      ITMAX = 25

      H = 0.01D0
      X = SC
      IT = 0
      TOL = 1.0D-9
  
  1   CONTINUE

      IT = IT + 1

      TMP = UC(N,C,K,0)

      CALL EQCSOLVE(N,C,K,X)
      VAL0 = X - UC(N,C,K,0)

      CALL EQCSOLVE(N,C,K,X+H)
      VAL1 = X+H - UC(N,C,K,0)
      CALL EQCSOLVE(N,C,K,X-H)
      VAL2 = X-H - UC(N,C,K,0)

      DIFF = (VAL2-VAL1)/(2.0D0*H)

      X = X + VAL0 / DIFF 

      IF(ABS(VAL0).GT.TOL.AND.IT.LT.ITMAX) GOTO 1 

      SC = X

      RETURN
      END

      SUBROUTINE EQITER2(N,C,K,SC)
C Find marginal social cost SC such that mass
C N+1 (late arrival) is just to become active
C Newton algorithm is used for this one dimensional problem. 
C The jacobian is approximated
C using finite differences.
      IMPLICIT NONE
      INTEGER N, C, LDA, IT,ITMAX
      PARAMETER(LDA=1000)
      DOUBLE PRECISION UC, K(N), H, VAL0,VAL1,VAL2,X,DIFF,SC,TOL
      DOUBLE PRECISION TMP

      ITMAX = 25

      H = 0.01D0
      X = SC
      IT = 0
      TOL = 1.0D-9
  
  1   CONTINUE

      IT = IT + 1

      K(N+1) = 0.0D0
      TMP = UC(N,C,K,N+1)

      CALL EQCSOLVE(N,C,K,X)
      VAL0 = X - UC(N,C,K,N+1)

      CALL EQCSOLVE(N,C,K,X+H)
      VAL1 = X+H - UC(N,C,K,N+1)
      CALL EQCSOLVE(N,C,K,X-H)
      VAL2 = X-H - UC(N,C,K,N+1)

      DIFF = (VAL2-VAL1)/(2.0D0*H)

      X = X + VAL0 / DIFF 

C     PRINT *, X, VAL0, DIFF, IT, K(1:3)
      IF(ABS(VAL0).GT.TOL.AND.IT.LT.ITMAX) GOTO 1 

      SC = X

      RETURN
      END

