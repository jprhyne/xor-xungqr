      ! Since we have a larft_ut with and without computing an inverse, we will
      ! time constructing T and applying to a matrix C (larfb). We do NOT test
      ! for correctness. This will be done in either the test directory or in
      ! the main lapack test suite.
      PROGRAM TIME_DLARFT
         ! Local Scalars
         INTEGER M,N,K,INFO,LWORK,NEG_ONE
         REAL START, FINISH
         ! Local Arrays
         DOUBLE PRECISION, ALLOCATABLE :: A(:,:), T(:,:), TAU(:),
     $      WORK(:),C(:,:)
         ! QR{TIMES,PERF} = [ OPT, REF_LVL2, REF_UT, UT_INV, UT_SOLVE ]
         DOUBLE PRECISION, DIMENSION(5) :: QRTIMES, QRPERF
         ! Intrinsic Functions
         INTRINSIC CPU_TIME,MAX
         ! External Functions
         ! External Subroutines
         ! Beginning of executable statements
         ! First, we read in the dimensions
!        WRITE(*,*)  "Provide N >= K for testing DLARF{B,T}",
!    $               "subroutines"
         READ(*,*) N, K
         NEG_ONE = -1
         ! We don't check this, but we will really only be checking
         ! for the cases of
         ! 1) N = K
         ! 2) N > K
         ! Allocate our memory
         ALLOCATE(A(N,N))
         ALLOCATE(C(N,N))
         ALLOCATE(T(K,K))
         ALLOCATE(TAU(N))
         ALLOCATE(WORK(1))
         ! Fill all our arrays with random data
         CALL RANDOM_NUMBER(A)
         CALL RANDOM_NUMBER(C)
         CALL RANDOM_NUMBER(T)
         CALL RANDOM_NUMBER(TAU)
         ! Perform our workspace queries.
         CALL DGEQRF(N, K, A, N, TAU, WORK, NEG_ONE, INFO)
         LWORK = INT(WORK(1))
         CALL DGELQF(K, N, A, N, TAU, WORK, NEG_ONE, INFO)
         LWORK = MAX(LWORK,INT(WORK(1)))
         DEALLOCATE(WORK)
         ALLOCATE(WORK(LWORK))
         ! Now, we factorize A
         ! QR
         CALL DGEQRF(N, K, A, N, TAU, WORK, LWORK, INFO)
         ! Now, we compute T then apply it to C
         ! AOCL
         CALL CPU_TIME(START)

         CALL DLARFT('F', 'C', N, K, A, N, TAU, T, K)

         CALL CPU_TIME(FINISH)
         QRTIMES(1) = FINISH - START
         ! Level 2 BLAS Implementation
         CALL CPU_TIME(START)

         CALL DLARFT_LVL2('F', 'C', N, K, A, N, TAU, T, K)

         CALL CPU_TIME(FINISH)
         QRTIMES(2) = FINISH - START
         ! Recursive (reference)
         CALL CPU_TIME(START)

         CALL DLARFT_REC('F', 'C', N, K, A, N, TAU, T, K)

         CALL CPU_TIME(FINISH)
         QRTIMES(3) = FINISH - START
         ! UT_INV
         CALL CPU_TIME(START)

         CALL DLARFT_UT('F', 'C', '1', N, K, A, N, TAU, T, K)

         CALL CPU_TIME(FINISH)
         QRTIMES(4) = FINISH - START
         ! UT_SOLVE
         CALL CPU_TIME(START)

         CALL DLARFT_UT('F', 'C', '2', N, K, A, N, TAU, T, K)

         CALL CPU_TIME(FINISH)
         QRTIMES(5) = FINISH - START
         ! Convert the execution times into performance metrics
         CALL QX_PERF(QRTIMES, QRPERF, N, K)
         WRITE(*,*) "QR TIMES:", QRTIMES
         WRITE(*,*) "QR PERF :", QRPERF
         !WRITE(*,*) QRTIMES
         !WRITE(*,*) QRPERF
         CALL DGELQF(K, N, A, N, TAU, WORK, LWORK, INFO)
         ! LQ
         CALL CPU_TIME(START)

         CALL DLARFT('F', 'R', N, K, A, N, TAU, T, K)

         CALL CPU_TIME(FINISH)
         QRTIMES(1) = FINISH - START
         ! Level 2 BLAS Implementation
         CALL CPU_TIME(START)

         CALL DLARFT_LVL2('F', 'R', N, K, A, N, TAU, T, K)

         CALL CPU_TIME(FINISH)
         QRTIMES(2) = FINISH - START
         ! Recursive (reference)
         CALL CPU_TIME(START)

         CALL DLARFT_REC('F', 'R', N, K, A, N, TAU, T, K)

         CALL CPU_TIME(FINISH)
         QRTIMES(3) = FINISH - START
         ! UT_INV
         CALL CPU_TIME(START)

         CALL DLARFT_UT('F', 'R', '1', N, K, A, N, TAU, T, K)

         CALL CPU_TIME(FINISH)
         QRTIMES(4) = FINISH - START
         ! UT_SOLVE
         CALL CPU_TIME(START)

         CALL DLARFT_UT('F', 'R', '2', N, K, A, N, TAU, T, K)

         CALL CPU_TIME(FINISH)
         QRTIMES(5) = FINISH - START
         ! Convert the execution times into performance metrics
         CALL QX_PERF(QRTIMES, QRPERF, N, K)
         WRITE(*,*) "LQ TIMES:", QRTIMES
         WRITE(*,*) "LQ PERF :", QRPERF
      END PROGRAM

      SUBROUTINE QX_PERF(QXTIMES, QXPERF, N, K)
         DOUBLE PRECISION, DIMENSION(5) :: QXTIMES, QXPERF
         INTEGER N, K
         DOUBLE PRECISION TMP,TIME,LARFB_FLOP, ZERO
         PARAMETER(TIME=1.0D+9, ZERO = 0.0D+0)
         ! Functions
         DOUBLE PRECISION LARFB0C2_PERF, LARFT_PERF, LARFT_UT_PERF
         EXTERNAL LARFB0C2_PERF, LARFT_PERF, LARFT_UT_PERF
         ! Compute performance metric for optimized and reference
         ! Store the flops for calls to larfb0c2 as these are the same
         LARFB_FLOP = ZERO !LARFB0C2_PERF(N, N, K, 1)
         TMP = LARFT_PERF(N, K) + LARFB_FLOP
         QXPERF(1:3) = TMP / (QXTIMES(1:3)*TIME)
         ! Compute performance metric for larft_ut with trtri call
         TMP = LARFT_UT_PERF(N, K, 1) + LARFB_FLOP
         QXPERF(4) = TMP / (QXTIMES(3)*TIME)
         ! Compute performance metric for larft_ut without trtri call
         TMP = LARFT_UT_PERF(N, K, 2) + LARFB_FLOP
         QXPERF(5) = TMP / (QXTIMES(4)*TIME)
      END SUBROUTINE

      DOUBLE PRECISION FUNCTION LARFT_PERF(N, K)
         DOUBLE PRECISION TWO, THREE, SIX
         INTEGER N, K

         PARAMETER(TWO=2.0D+0, THREE=3.0D+0, SIX=6.0D+0)

         LARFT_PERF = SIX*N*K*K - SIX*N*K - TWO*K*K*K + THREE*K*K - K
         LARFT_PERF = LARFT_PERF / SIX
      END FUNCTION

      DOUBLE PRECISION FUNCTION LARFT_UT_PERF(N, K, VAR)
         DOUBLE PRECISION TWO, THREE, FOUR, SIX, TMP
         INTEGER N, K, VAR

         PARAMETER(TWO=2.0D+0, THREE=3.0D+0, FOUR=4.0D+0, SIX=6.0D+0)
         ! var.eq.1 means we call trtri
         ! var.eq.2 means we don't call trtri
         IF (VAR.EQ.1) THEN
            TMP = SIX*N*K*K + SIX*N*K - TWO*K*K*K - THREE*K*K + TWO*K
         ELSE IF (VAR.EQ.2) THEN
            TMP = SIX*N*K*K + SIX*N*K - FOUR*K*K*K - THREE*K*K - TWO*K
         END IF
         LARFT_UT_PERF = TMP / SIX
      END FUNCTION

      DOUBLE PRECISION FUNCTION LARFB0C2_PERF(M, N, K, VAR)
         DOUBLE PRECISION TWO, FOUR
         INTEGER M, N, K, VAR

         PARAMETER(TWO=2.0D+0, FOUR=4.0D+0)

         IF (VAR.EQ.1) THEN ! QX
            LARFB0C2_PERF=FOUR*M*N*K - TWO*N*K*K
         ELSE IF (VAR.EQ.2) THEN ! XQ
            LARFB0C2_PERF=FOUR*M*N*K - TWO*M*K*K
         END IF
      END FUNCTION
