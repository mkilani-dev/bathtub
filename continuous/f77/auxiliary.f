
      SUBROUTINE READPARAM(FILENAME)
      IMPLICIT NONE
C Read file for parameter values and store these in
C common blocks

      INTEGER N, LENGTH, IND, NGRID 
      DOUBLE PRECISION R0,R1,S0,S1,VF,K,L,P,TMAX
      DOUBLE PRECISION M0, MBAR
      INTEGER          ITMAX
      CHARACTER*10 FILENAME
      CHARACTER*20 LINE
      CHARACTER*6 VARNAME
      DOUBLE PRECISION VARVALUE

      COMMON /MODELPARAM/  R0,R1,S0,S1,VF,K,L,P,TMAX
      COMMON /OTHERPARAM/  M0, MBAR, ITMAX

      OPEN(UNIT=1,FILE=FILENAME)
 
 1    READ(1,'(a20)',END=3) LINE
      LENGTH=SCAN(LINE,'#',.FALSE.)
      IF(LENGTH.GT.0) GOTO 1
      LENGTH=SCAN(LINE,'=',.FALSE.)
      IF(LENGTH.EQ.0) GOTO 1
      READ(LINE(1:LENGTH-1),'(A6)') VARNAME
      READ(LINE(LENGTH+1:40),*) VARVALUE
 
      IF(VARNAME.EQ.'VF') THEN
         VF=VARVALUE
      ELSEIF(VARNAME.EQ.'P') THEN
         P=VARVALUE
      ELSEIF(VARNAME.EQ.'L') THEN
         L=VARVALUE
      ELSEIF(VARNAME.EQ.'K') THEN
         K=VARVALUE
      ELSEIF(VARNAME.EQ.'TMAX') THEN
         TMAX=VARVALUE
      ELSEIF(VARNAME.EQ.'R0') THEN
         R0=VARVALUE
      ELSEIF(VARNAME.EQ.'R1') THEN
         R1=VARVALUE
      ELSEIF(VARNAME.EQ.'S0') THEN
         S0=VARVALUE
      ELSEIF(VARNAME.EQ.'S1') THEN
         S1=VARVALUE 
      ELSEIF(VARNAME.EQ.'M0') THEN
         M0=VARVALUE 
      ELSEIF(VARNAME.EQ.'MBAR') THEN
         MBAR=VARVALUE 
      ELSEIF(VARNAME.EQ.'ITMAX') THEN
         ITMAX=NINT(VARVALUE)
      ENDIF

      GOTO 1

 3    CONTINUE
      CLOSE(1)

      RETURN
      END

C ==================================================
C Descriptive values for given vector of entry rates
C ==================================================
      SUBROUTINE STATDES(N, E, DAT, IDAT)
      IMPLICIT NONE
      INTEGER          I, N, NL, NMAX, NE
      PARAMETER        (NMAX=10000)
      DOUBLE PRECISION U, E(N), EE(N)
      DOUBLE PRECISION DAT(NMAX)
      INTEGER          IDAT(NMAX)
C Variables used for intermediate computations
      DOUBLE PRECISION V(NMAX), T(NMAX), TT(NMAX)
      DOUBLE PRECISION UE(NMAX), UX(NMAX)
      DOUBLE PRECISION SUME(NMAX),CUME(NMAX)
C Parameters of the problem
      DOUBLE PRECISION R0, R1, S0, S1, VF, K, L, P, TMAX 
      COMMON /MODELPARAM/  R0, R1, S0, S1, VF, K, L, P, TMAX

      NL = IDAT(2) 
      NE = IDAT(3)


C Compute travel speed, travel time and clock-time on each interval
      I=1
      TT(I) = DAT(I) / VF
      SUME(I)=E(I)
      V(I) = VF * ( 1.0D0 - P * SUME(I) / K )
      T(I) = ( DAT(I+1) - DAT(I) ) / V(I)
      TT(I+1) = TT(I) + T(I)
      DO 1, I = 2, N + NL - 1
        IF(I.LE.NL) THEN
          SUME(I) = SUME(I-1) + E(I)
        ELSEIF (I.LE.N) THEN
          SUME(I) = SUME(I-1) + E(I) - E(I-NL)
        ELSE
          SUME(I) = SUME(I-1)        - E(I-NL)
        ENDIF
        V(I) = VF * ( 1.0D0 - P * SUME(I) / K )
        T(I) = ( DAT(I+1) - DAT(I) ) / V(I)
        TT(I+1) = TT(I) + T(I)
 1    CONTINUE 

C Compute entry and exit sub-utilities
      DO 20, I = 1, N
         UE(I) = R0 * DLOG( R1 *          TT(I     ) ) * E(I) * P
         UX(I) = S0 * DLOG( S1 * ( TMAX - TT(I + NL))) * E(I) * P
 20   CONTINUE


      EE = 0.0D0
      CUME(1) = E(1)
      DO 30, I = 1, N+NL-1
      IF(MOD(I-1,10).EQ.0) PRINT'(A5,10A9)',"I","M","e","E","e-x",
     1     "t","v","Te","Ue","Tx","Ux"
      IF(I.LE.N-1) THEN
        CUME(I+1)=CUME(I)+E(I+1)
      ELSE
        CUME(I+1)=CUME(I)
      ENDIF
      IF(I.LE.N) EE(I) = E(I)
      PRINT'(I5,10F9.3)',I,DAT(I),EE(I),CUME(I),SUME(I),
     1                  T(I),V(I),TT(I),UE(I),TT(I+NL),UX(I)
 30   CONTINUE
      I=N+NL
      PRINT'(I5,10F9.3)',I,DAT(I),EE(I),CUME(I),0.0D0,
     1                  T(I),VF,TT(I),UE(I),TT(I+NL),UX(I) 


      RETURN
      END

      SUBROUTINE CHECK_GRAD(N, E, M, CON, DAT, IDAT)
      IMPLICIT NONE
C IDAT(2) = NL
C DAT (3) = K
C DAT (4) = P
C DAT (5) = TMAX 
      INTEGER          N, M, NMAX, NE, NL
      PARAMETER (NMAX=10000)
      DOUBLE PRECISION CON(M), E(N)
      DOUBLE PRECISION DAT(NMAX)
      INTEGER          IDAT(NMAX) 
      DOUBLE PRECISION D_T(N+IDAT(2)-1), D_TT(N+IDAT(2),N)
      DOUBLE PRECISION V(N+IDAT(2)-1), T(N+IDAT(2)-1), TT(N+IDAT(2))

      DOUBLE PRECISION SUME
      INTEGER          I, J

C Parameters of the problem
      DOUBLE PRECISION R0, R1, S0, S1, VF, K, L, P, TMAX 
      COMMON /MODELPARAM/  R0, R1, S0, S1, VF, K, L, P, TMAX


      NL = IDAT(2)
      NE = IDAT(3)

C Constraints 1 to NE - NL + 1

      DO 10, I = 1, NE - NL + 1
      SUME = 0.0D0
      DO 15, J = I, I + NL -1
      SUME = SUME + E(J)
15    CONTINUE
      CON(I) = SUME + E(NE+I) - K / P
10    CONTINUE

C Constraint NE - NL + 1 + 1 (copied from UTIL, except last line) 
      I=1
      TT(I) = DAT(I) / VF
      SUME=E(I)
      V(I) = VF * ( 1.0D0 - P * SUME / K )
      T(I) = ( DAT(I+1) - DAT(I) ) / V(I)
      TT(I+1) = TT(I) + T(I)
      DO 20, I = 2, N + NL - 1
        IF(I.LE.NL) THEN
          SUME = SUME + E(I)
        ELSEIF (I.LE.NE) THEN
          SUME = SUME + E(I) - E(I-NL)
        ELSE
          SUME = SUME        - E(I-NL)
        ENDIF
        V(I) = VF * ( 1.0D0 - P * SUME / K )
        T(I) = ( DAT(I+1) - DAT(I) ) / V(I)
        TT(I+1) = TT(I) + T(I)
 20   CONTINUE 
      CON(NE-NL+1+1) = TT(NE+NL) + E(NE+NE-NL+1+1) - TMAX


C Constraint NE + NE - NL + 1 + 1 + 1
      SUME = 0.0D0
      DO 30, I = 1, NE
      SUME = SUME + E(I)
 30   CONTINUE 
      CON(NE-NL+3) = SUME - 1.0D0 

      RETURN
      END
