
      PROGRAM MASS
C This program computes optimum mass entries in the bathtub model
C 
C The main program and subroutine are in this file. Routines from
C reading dat input and displaying output are in put in file
C "auxiliary.f".
C 
C The main program call a first subroutines to read data input and then
C calls the subroutine that computes the optimum entry rates.
C
C A SUMMARY OF THE SUBROUTINES 
C
C UTIL      : computes the objective function
C D_UTIL    : computes the gradient of the objective function
C CONS      : computes the value of the constraint. This version 
C             has one constraint (the sum of the entry rates is equal to
C             one
C D_CONS    : computes the gradient of the constraint
C SPLITLINE : for input, it takes the current approximation line with
C             respect variable m and breaks and evenly subdivides the largest
C             intervals into two intervals. Bu doing so we intend to get
C             a smoother approximation of the solution
C
C For the optimization step we use the fortran version IPOPT (the large-scale 
C optimization package): 
C
C https://www.coin-or.org/Ipopt/ipopt-fortran.html
C 
C TO RUN THIS PROGRAM

C - install IPOPT (it requires LAPACK and BLAS)
C - compile through a command like 
C       gfortran -o massn massn.f auxiliary.f -llapack -lblas -lipopt  
C   by providing the correct path to the libraries
C - set data in file "params.dat"
C - run the program
C       ./massn
C
C The process is little tedious, in particular for IPOPT. Another
C alternative, is that I can provide you with an access to my personal
C server and you connect through SSH.
C
C
C KNOWN ISSUES AND POSSIBLE IMPROVEMENTS
C
C 1. This version does not explicitly set constraint on :
C    - last arrival time should be smaller than TMAX
C    - traffic density should be smaller than road jam capacity
C    and relies on the optimization to carry on implicitly these
C    constraints. This does not occur in all cases and fine tuning of
C    parameter M0, MBAR and ITMAX is frequently required for convergence.
C
C 2. M0 and MBAR are fixed. It is suggested to provide a large window in
C    the start. The solution will entail a first phase where without
C    entries. Similarly it will entail a last phase where entries are
C    equal to zero. A new of the computation can then be considered by
C    zooming a new window with higher M0 and smaller MBAR. This process
C    can be repeated until satisfactory values are obtained.
C
C 3. The Hessian is not provided to the optimizer. It uses finite
C    difference to approximate it (IQUASI = 6, in the options passed to
C    IPOPT).
C
C 4. The program does not produce any graphical output. For instance,
C    output table can be exported to a text file (redirect the output
C    using the unix standard "./massn > outfile") and then use gnuplot
C    or another program to make the plots (after deleting useless
C    lines).
C
C These issues will be addressed in the updates
C
C
C NOTATION THROUGH AN EXAMPLE   
C  - number of entry (and exit) points: N = 10
C  - number of entry points before within a cycle:  NL = 2
C  - total number of breakpoints: N + NL = 12
C
C (0)    t1     t2    t3     t4    t5     t6     t7    t8     t9   t10    t11      
C (1)  P----S-------P----S-------P----S-------P-----S-------P----S-----P--------S--- 
C
C (2)  0    Mb-4L   L    Mb-3L   2L   Mb-2L   3L    Mb-L    4L   Mb    (4+1)L   Mb+L   
C
C (3)  1    2       3    4       5    6       7     8       9    10    11       12 
C
C (4)  e1   e2      e3   e4      e5   e6      e7    e8      e9   e10   
C
C (5)               e1   e2      e3   e4      e5    e6      e7   e8    e9       e10   
C 
C (6)  T1   T2      T3   T4      T5   T6      T7    T8      T9   T10   T11      T12
C
C
C  (0) t_i clock-time from m_i to m_(i+1)
C  (1) Primary (P) and secondary (S) breakpoints
C  (3) indices for variable m
C  (4) entries at each breakpoint m_i (at each breakpoint)
C  (5) exits at each breakpoint m_i (at each breakpoint)
C  (6) T_i clock-time at mile m_i 

      IMPLICIT NONE
      
      DOUBLE PRECISION UT
      DOUBLE PRECISION R0, R1, S0, S1, VF, K, L, P, TMAX 
      DOUBLE PRECISION M0, MBAR
      INTEGER          ITMAX
      COMMON /MODELPARAM/  R0, R1, S0, S1, VF, K, L, P, TMAX 
      COMMON /OTHERPARAM/  M0, MBAR, ITMAX

C Read parameters from file "params.dat"
      CALL READPARAM("params.dat")
      PRINT '(11F9.2)', R0,R1,S0,S1,VF,K,L,P,TMAX,M0,MBAR
C Call subroutine UTILOPT to compute the optimum entry rates
      CALL UTILOPT ( M0, MBAR , UT, ITMAX)

C Print the value of the objective function
      PRINT*, UT

      END


      SUBROUTINE UTILOPT (M0,MBAR,F,ITMAX)
C Computes entry masses at points m_i that maximize aggregate
C utility given that M0 and MBAR are fixed
      IMPLICIT NONE

      DOUBLE PRECISION M0, MBAR
      INTEGER          I, NCYC, NL, NMAX, IT, ITMAX
      PARAMETER        (NMAX=10000)
      DOUBLE PRECISION U

C Parameters of the problem
      DOUBLE PRECISION R0, R1, S0, S1, VF, K, L, P, TMAX 

C Variable definition for IPOPT
      INTEGER          LRW,  LIW
      PARAMETER        ( LRW = 10000, LIW = 10000)
      INTEGER          IW(LIW)
      DOUBLE PRECISION RW(LRW)
      DOUBLE PRECISION IPOPT
      INTEGER          N, M
      DOUBLE PRECISION LAM(NMAX)
      DOUBLE PRECISION C(NMAX)
      DOUBLE PRECISION G(NMAX)
      DOUBLE PRECISION E(NMAX)
      INTEGER          NLB, NUB
      INTEGER          ILB(NMAX), IUB(NMAX)
      DOUBLE PRECISION BNDSL(NMAX), BNDSU(NMAX)
      DOUBLE PRECISION V_L(NMAX), V_U(NMAX)
      DOUBLE PRECISION DAT  (NMAX)
      INTEGER          IDAT (NMAX)
      INTEGER          NARGS
      DOUBLE PRECISION ARGS (50)
      CHARACTER*20     CARGS(50)
      INTEGER          IERR, ITER
      DOUBLE PRECISION F

      EXTERNAL       UTIL, CONS, D_UTIL, D_CONS
      EXTERNAL       EV_H_DUMMY,EV_HLV_DUMMY, EV_HOV_DUMMY, EV_HCV_DUMMY
C End of definition of IPOPT variables

      COMMON /MODELPARAM/  R0, R1, S0, S1, VF, K, L, P, TMAX
 
      NCYC= FLOOR ( (MBAR-M0) / L )  
      N = 2 * (  NCYC+ 1 )
      IF(MBAR.GT.L) THEN
        NL = 2
      ELSE
        NL =1
      ENDIF
      IDAT(1) = NCYC
      IDAT(2) = NL

C Compute the breakpoints. This is the first approximation that will be
C used.
      DO 10, I = 1,  NCYC+ NL
C Primary breakpoints
        DAT(2*I-1) = M0   + DBLE( I - 1 ) * L
C Secondary breakpoints. Notice that : 
C       MBAR + L -  NCYCL + (I-1) * L = MBAR - (NCYC-I) * L 
        DAT(2*I)   = MBAR - DBLE (  NCYC- I + 1 ) * L
 10   CONTINUE

C Compute a simple initial solution
      DO 20, I = 1, N
        E(I) = ( DAT(I+1) - DAT(I) ) / ( (  NCYC + 1 ) * L ) 
 20   CONTINUE

C Set algorithmic parameters (for IPOPT)
       NARGS   = 1
       ARGS(1) = 1.D-8
       CARGS   = 'dtol' 


C MAIN LOOP
C It starts by optimizing over the initial breakpoints, then 
C uses the obtained solution to increase the number of approximation
C points. Variable ITMAX indicates the number of iterations (how many
C time the approximation line) is subdivided.

      IT = 1
 100  CONTINUE

C Number of variable with lower bounds
      NLB = N 
C Indices of variables with lower bounds
      DO 1, I = 1, N
      ILB(I) = I
 1    CONTINUE
C Values of lower bounds
      DO 2, I = 1, N
      BNDSL(I) = 1.0D-6
 2    CONTINUE
C Number of variables with upper bounds
      NUB = 0
C i.e. no uppers bounds on the entry rates

C There is only one constraint on the problem (the sum of all entry
C rates is one)
      M = 1
      PRINT* , N, NL, IT

C  Call IPOPT (the optimizer)
      F = - IPOPT(N, E, M, NLB, ILB, BNDSL, NUB, IUB, BNDSU, V_L, V_U,
     1   LAM, C, LRW, RW, LIW, IW, ITER, IERR, UTIL, CONS, D_UTIL,
     1   D_CONS, EV_H_DUMMY, EV_HLV_DUMMY, EV_HOV_DUMMY, EV_HCV_DUMMY,
     1   DAT, IDAT, NARGS, ARGS, CARGS)

      PRINT '(A20,F9.3)' , 'Value of objective function : ', F 

       PRINT '(F9.5)', SUM(E) 

      IT = IT + 1 

      IF(IT.LE.ITMAX) THEN
        CALL SPLITLINE(N, DAT,NL,E,0.80D0)
        IDAT(2) = NL
        GOTO 100
      ENDIF

C Print the last solution obtained
      CALL STATDES(N, E, DAT, IDAT)

      RETURN
      END

C ======================================
C Computation of the objective function
C (aggregate utility)
C ======================================
      SUBROUTINE UTIL(N, E, U, DAT, IDAT)
      IMPLICIT NONE
      INTEGER          I, N, NL, NMAX
      PARAMETER        (NMAX=10000)
      DOUBLE PRECISION U, E(N)
      DOUBLE PRECISION DAT(NMAX)
      INTEGER          IDAT(NMAX)
C Variables used for intermediate computations
      DOUBLE PRECISION V(N+IDAT(2)-1), T(N+IDAT(2)-1), TT(N+IDAT(2))
      DOUBLE PRECISION UE, UX
      DOUBLE PRECISION SUME
C Parameters of the problem
      DOUBLE PRECISION R0, R1, S0, S1, VF, K, L, P, TMAX 
      COMMON /MODELPARAM/  R0, R1, S0, S1, VF, K, L, P, TMAX

      NL = IDAT(2) 

C Compute travel speed, travel time and clock-time on each interval
      I=1
      TT(I) = DAT(I) / VF
      SUME=E(I)
      V(I) = VF * ( 1.0D0 - P * SUME / K )
      T(I) = ( DAT(I+1) - DAT(I) ) / V(I)
C     PRINT*, V(I), VF, P, K, SUME,T(I),TT(I),DAT(2)- DAT(1)
      TT(I+1) = TT(I) + T(I)
      DO 1, I = 2, N + NL - 1
        IF(I.LE.NL) THEN
          SUME = SUME + E(I)
        ELSEIF (I.LE.N) THEN
          SUME = SUME + E(I) - E(I-NL)
        ELSE
          SUME = SUME        - E(I-NL)
        ENDIF
        V(I) = MAX(0.D0,VF * ( 1.0D0 - P * SUME / K ))
        T(I) = ( DAT(I+1) - DAT(I) ) / V(I)
        TT(I+1) = TT(I) + T(I)
C       PRINT '(I5,5F9.3)', I,DAT(I),V(I), T(I), TT(I+1), SUME
 1    CONTINUE 

C Compute entry and exit sub-utilities
      UE = 0.0D0 
      UX = 0.0D0 
      DO 20, I = 1, N
         UE = UE + R0 * DLOG( R1 *          TT(I     ) ) * E(I) * P
         UX = UX + S0 * DLOG( S1 * ( TMAX - TT(I + NL))) * E(I) * P
C     PRINT '(2I5,4F9.3)', I, NL,TMAX,TT(I+NL),UE, UX
 20   CONTINUE


C AGREGATE UTILITY 
C Notice that we consider the opposite of the utility since we set the
C problem as a minimization program
      U = - ( UE + UX )  


      RETURN
      END


C =====================================================
C Computation of the gradient of the objective function
C =====================================================
      SUBROUTINE D_UTIL( N, E, G, DAT, IDAT)
      IMPLICIT NONE

      INTEGER          N, NL, NMAX
      PARAMETER        (NMAX=10000)
      DOUBLE PRECISION G(N), E(N)
      DOUBLE PRECISION DAT(NMAX)
      INTEGER          IDAT(NMAX)
      DOUBLE PRECISION D_T(N+IDAT(2)-1), D_TT(N+IDAT(2),N)
      DOUBLE PRECISION D_UE(N,N), D_UX(N,N)
      DOUBLE PRECISION V(N+IDAT(2)-1), T(N+IDAT(2)-1), TT(N+IDAT(2))
      DOUBLE PRECISION SUME
      INTEGER          I, J, II
C Parameters of the problem
      DOUBLE PRECISION R0, R1, S0, S1, VF, K, L, P, TMAX 
      COMMON /MODELPARAM/  R0, R1, S0, S1, VF, K, L, P, TMAX
      
      NL = IDAT(2)
C Compute travel speed on each interval
C This step is replicated from UTIL subroutine.
      I = 1
      TT(I)  = DAT(I) / VF
      SUME   = E(I)
      V(I)   = MAX(0.0D0,VF * ( 1.D0 - P * SUME / K ))
      T(I)   = ( DAT(I+1) - DAT(I) ) / V(I)
      TT(I+1) = TT(I) + T(I)
C     PRINT*, V(I), VF, P, K, SUME,T(I),TT(I),DAT(2)- DAT(1)
      DO 1, I = 2, N + NL - 1
        IF(I.LE.NL) THEN
          SUME = SUME + E(I)
        ELSEIF (I.LE.N) THEN
          SUME = SUME + E(I) - E(I-NL)
        ELSE
          SUME = SUME        - E(I-NL)
        ENDIF
        V(I) = VF * ( 1.D0 - P * SUME / K )
        T(I) = ( DAT(I+1) - DAT(I) ) / V(I)
        TT(I+1) = TT(I) + T(I)
C     PRINT*, V(I), VF, P, K, SUME,T(I),TT(I),DAT(I+1), DAT(I)
 1    CONTINUE 

C Compute the derivative
C
C   d t_i    m_(i+1) - m_i     vf * P
C   ----- =  -------------  *  ------
C   d e_i        v_i^2           K
C
      DO 10, I = 1, N + NL - 1
      D_T(I) = (DAT(I+1) - DAT(I))/V(I)**2 * VF*P/K
 10   CONTINUE 
      
C Compute the derivative
C
C                  /
C                 |   0                            if  j <= i
C                 |      
C                 |
C      d T_j      |  d t_i           d t_(j-1)
C      -----  =  /   -----  + ... +  ---------     if  i < j <=i+NL 
C      d e_i     \   d e_i           d e_(j-1)
C                 |      
C                 |
C                 |  d T_(i+NL)    
C                 |  ----------                    if  j > i + NL
C                 |   d e_i         
C                  \
C 
      
      D_TT = 0.0D0
      DO 20, I = 1, N
      DO 20, J = I+1, N+NL
        IF (J.LE.I) THEN
          D_TT(J,I) = 0.0D0
        ELSEIF (J.LE.I+NL) THEN
          DO 22, II = I, J-1
          D_TT(J,I) = D_TT(J,I) + D_T(II)
22        CONTINUE
        ELSE
          D_TT(J,I) = D_TT(I+NL,I)
        ENDIF 
20    CONTINUE 

C IMPACTS ON ENTRY SUBUTILITIES
C An increase of the entry rate at T_i increases the size of the group
C itself and delays the entry time for all groups that enter at T_j with
C j > i.
C Users who enter at T_N do not delay the entry of any other group. The
C structure of the loop automatically takes into account this case since.
C Indeed, for I=N, J varies varies from N+1 to N and, then, the loop is
C not executed.
      D_UE = 0.0D0
      DO 30, I =   1, N 
      D_UE(I, I) = R0 * DLOG( R1 * TT(I) ) * P
      DO 30, J = I+1, N
      D_UE(J, I) = R0 * D_TT(J,I) / TT(J) * E(J) * P 
30    CONTINUE 

C IMPACTS ON EXIT SUBUTILITIES
C An increase of the size of the group of users entering at T_i
C has two impacts: (1) it delays the arrival of all groups entering at 
C T_{max(1,i-(NL-1))}, including the group itself, and (2) it increases
C the size of the group of those entering at T_i. 
      D_UX = 0.0D0 
      DO 40, I  = 1, N 
      DO 40, J  = MAX(I-(NL-1), 1), N
      D_UX(J, I) = - S0 * D_TT(J+NL,I) / (TMAX - TT(J+NL)) * E(J) * P 
40    CONTINUE
      DO 45, I  = 1, N
      D_UX(I, I) = D_UX(I, I) + S0 * DLOG( S1 * (TMAX - TT(I+NL)) ) * P
45    CONTINUE 

C TOTAL IMPACTS (the gradient)
C Notice that we compute the gradient of - UTIL since we set
C the problem as a minimization program.
      G = 0.0D0
      DO 50, I = 1, N
      DO 50, J = 1, N 
      G(I) = G(I) - D_UE(J,I) - D_UX(J,I)
C     PRINT '(7F9.3)', G(I), E(1:6) 
 50   CONTINUE 

      RETURN
      END 


C =========================================
C Computation of the equality constraints
C In this problem, there is only one constraint
C stating that the sum of all entry rates 
C is equal to one :
C
C      e_1 + e_2 + .. + e_n = 1 
C
C ========================================= 
      SUBROUTINE CONS(N, E, M, C, DAT, IDAT)
      IMPLICIT NONE

      INTEGER          N, M, NMAX
      PARAMETER (NMAX=10000)
      DOUBLE PRECISION C(M), E(N)
      DOUBLE PRECISION DAT(NMAX)
      INTEGER          IDAT(NMAX)


      DOUBLE PRECISION SUME
      INTEGER          I

      SUME = 0.D0
      DO 10, I = 1, N
        SUME = SUME + E(I)
 10   CONTINUE

      C(1) = SUME - 1.D0

      RETURN
      END
 

C =============================================
C Computation of the Jacobian of the constraint
C =============================================
      SUBROUTINE D_CONS(TASK, N, E, NZ, A, ACON, AVAR, DAT, IDAT)
      IMPLICIT NONE
      INTEGER          TASK, N, NZ
      INTEGER          NMAX, I
      PARAMETER (NMAX=10000)
      DOUBLE PRECISION E(N), A(NZ)
      INTEGER          ACON(NZ), AVAR(NZ)
      DOUBLE PRECISION DAT(NMAX)
      INTEGER          IDAT(NMAX)
      
      IF(TASK.EQ.0) THEN
        NZ = N
      ELSE
        DO 10, I = 1, N
        AVAR(I) = I
        ACON(I) = 1
10      CONTINUE 
        DO 20, I = 1, N
        A(I) = 1.D0
20      CONTINUE
      ENDIF

      RETURN
      END


C =============================================
C Split interval
C =============================================
      SUBROUTINE SPLITLINE(N, M, NL, E, THR)
      IMPLICIT NONE
C Increases the number of approximation points by breaking the largest
C intervals.
C N   the number of entry rates, or mile-points where entry occurs
C M   table containing the approximation points
C NL  index the last element between0 and L: we have M(NL+1) = L  
C THR intervals with length larger than THR * X are swplit into two
C     intervals, where X is the the length of the largest interval

      INTEGER          NL, N, I, J, NMAX
      PARAMETER        (NMAX=10000) 
      DOUBLE PRECISION M(NMAX), E(NMAX)
      DOUBLE PRECISION THR
      DOUBLE PRECISION M2(NMAX), E2(NMAX)
      INTEGER          N2,NPTS,NL2,NTOT
      DOUBLE PRECISION Y(NMAX), YMAX

      IF(N.LT.2) THEN
        PRINT*, "A table of legth 2 at least is expected"
        STOP
      ENDIF
      IF(THR.LE.0.0D0.OR.THR.GE.1.0D0) THEN
        PRINT*, "THR should be between 1/2 and 1."
        STOP
      ENDIF 


C Compute the length of the largest interval
      YMAX = M(2) - M(1)
      DO 10, I = 1, NL
        Y(I) = M(I+1) - M(I)
        IF(Y(I).GT.YMAX) YMAX = Y(I)
 10   CONTINUE 

C Store some values
      N2 = N
      NL2 = NL
      NTOT = N + NL
      E2 = 0.0D0

C Split intervals with length higher than
C THR * YMAX and update indexes
      J = 0
      DO 20, I = 1, N+NL-1 
        IF((M(I+1)-M(I)).GT.(THR*YMAX)) THEN
          M2(I+J  ) = M(I)
          M2(I+J+1) = M(I) + (M(I+1) - M(I)) / 2.0D0
          J  = J + 1
          IF(I.LE.NL) THEN
            E2(I+J  ) = E(I)
            E2(I+J-1) = 0.0D0
            NL2  = NL2 + 1
            N2   = N2  + 1
            NTOT = NTOT + 1
          ELSEIF(I.LT.N) THEN
            E2(I+J  ) = E(I)
            E2(I+J-1) = 0.0D0
            N2   = N2  + 1
            NTOT = NTOT + 1
          ELSE
            NTOT = NTOT + 1
          ENDIF
        ELSE
          M2(I+J) = M(I)
          E2(I+J) = E(I)
        ENDIF 
 20   CONTINUE 
      E2(N2  ) = E(N)
      M2(NTOT) = M(N+NL)

C Set valriables to new values
      M  = M2
      E  = E2
      N  = N2
      NL = NL2


      RETURN
      END

