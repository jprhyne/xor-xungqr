c Label definitions
c 10 is for the first set of timing statements (optimized larft)
c 20 is for the second set of timing statements (my larft)
      PROGRAM PROFILE_ORGQR
      INTEGER     M, N, K, NB, LWORK, INFO, NEG_ONE
      REAL START, FINISH
      DOUBLE PRECISION  ZERO

      DOUBLE PRECISION, ALLOCATABLE :: V(:,:), Q(:,:),
     $   WORK(:), TAU(:)

      ! Col 1 is with optimized routines
      ! Col 2 is with my recursive dlarft
      DOUBLE PRECISION, DIMENSION(4,2) :: TIMEMAT
      !

      INTRINSIC CPU_TIME

      NEG_ONE = -1
      ZERO = 0.0D+0

      WRITE(*,*) "Provide M,N,K for ORGQR"
      READ(*,*) M, N, K

      ALLOCATE(V(M,K))
      ALLOCATE(Q(M,N))
      ALLOCATE(TAU(K))
      ALLOCATE(WORK(1))

      CALL RANDOM_NUMBER(V)
      CALL RANDOM_NUMBER(TAU)

      CALL DGEQRF(M, K, V, M, TAU, WORK, NEG_ONE, INFO)
      LWORK = INT(WORK(1))
      DEALLOCATE(WORK)
      ALLOCATE(WORK(LWORK))

10    CALL DGEQRF(M, K, V, M, TAU, WORK, LWORK, INFO)
      CALL DLASET('All', M, N, ZERO, ZERO, Q, M)
      CALL DLACPY('All', M, K, V, M, Q, M)
      ! construct the triangular factor
      CALL CPU_TIME(START)

      CALL DLARFT('F','C', M, K, Q, M, TAU, Q, M)

      CALL CPU_TIME(FINISH)
      TIMEMAT(1,1) = FINISH - START
      ! Apply H to the trailing n-k columns of Q
      CALL CPU_TIME(START)

      CALL DLARFB0C2('I', 'B', 'L', 'N', 'F', 'C', M, N-K, K, Q, M,
     $      Q, M, Q(1,K+1), M)

      CALL CPU_TIME(FINISH)
      TIMEMAT(2,1) = FINISH - START
      ! Apply H to the first k columns of Q
      CALL CPU_TIME(START)

      CALL DORGKR('1', M,K,Q,M)

      CALL CPU_TIME(FINISH)
      TIMEMAT(3,1) = FINISH - START

      ! Set Q back to 0
      CALL DLASET('All', M, N, ZERO, ZERO, Q, M)
      ! Copy V back into Q
      CALL DLACPY('All', M, K, V, M, Q, M)

      CALL CPU_TIME(START)

      CALL DORGQR_OPT_K(M, N, K, Q, M, TAU, WORK, LWORK, INFO)

      CALL CPU_TIME(FINISH)
      TIMEMAT(4,1) = FINISH - START

20    CALL DGEQRF(M, K, V, M, TAU, WORK, LWORK, INFO)
      CALL DLASET('All', M, N, ZERO, ZERO, Q, M)
      CALL DLACPY('All', M, K, V, M, Q, M)
      ! construct the triangular factor
      CALL CPU_TIME(START)

      CALL DLARFT_REC('F','C', M, K, Q, M, TAU, Q, M)

      CALL CPU_TIME(FINISH)
      TIMEMAT(1,2) = FINISH - START
      ! Apply H to the trailing n-k columns of Q
      CALL CPU_TIME(START)

      CALL DLARFB0C2('I', 'B', 'L', 'N', 'F', 'C', M, N-K, K, Q, M,
     $      Q, M, Q(1,K+1), M)

      CALL CPU_TIME(FINISH)
      TIMEMAT(2,2) = FINISH - START
      ! Apply H to the first k columns of Q
      CALL CPU_TIME(START)

      CALL DORGKR('1', M,K,Q,M)

      CALL CPU_TIME(FINISH)
      TIMEMAT(3,2) = FINISH - START

      ! Set Q back to 0
      CALL DLASET('All', M, N, ZERO, ZERO, Q, M)
      ! Copy V back into Q
      CALL DLACPY('All', M, K, V, M, Q, M)

      CALL CPU_TIME(START)

      CALL DORGQR_REC_K(M, N, K, Q, M, TAU, WORK, LWORK, INFO)

      CALL CPU_TIME(FINISH)
      TIMEMAT(4,2) = FINISH - START

      WRITE(*,*) TIMEMAT(1:4,1)
      WRITE(*,*) TIMEMAT(1:4,2)
      END PROGRAM
