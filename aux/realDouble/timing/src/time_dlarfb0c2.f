      PROGRAM TIME_LARFB
         ! Local scalars
         INTEGER M, N, K, MK, INFO, NEG_ONE, LWORK
         REAL START, FINISH
         ! Local Arrays
         DOUBLE PRECISION, ALLOCATABLE :: V(:,:), A(:,:), T(:,:),
     $      C(:,:), TAU(:), WORK(:), T_INV(:,:), WORK_MAT(:,:)

         !{TIME,PERF}VEC = [ AOCL (LARFB), REFERENCE (LARFB),
         !                   LARFB0C2_MULT LARFB0C2_SOLVE ]
         DOUBLE PRECISION, DIMENSION(4) :: TIMEVEC, PERFVEC

         INTRINSIC   CPU_TIME

         NEG_ONE = -1

         !WRITE(*,*) "Provide M,N,K for LARFB{,0C2}"
         READ(*,*) M, N, K

         ALLOCATE(A(M,K))
         ALLOCATE(V(M,K))
         ALLOCATE(T(K,K))
         ALLOCATE(T_INV(K,K))
         ALLOCATE(C(M,N))
         ALLOCATE(TAU(M))
         ALLOCATE(WORK(1))
         ALLOCATE(WORK_MAT(M,K))

         CALL RANDOM_NUMBER(A)
         CALL RANDOM_NUMBER(V)
         CALL RANDOM_NUMBER(T)
         CALL RANDOM_NUMBER(C)
         CALL RANDOM_NUMBER(TAU)

         ! Copy A into V
         CALL DLACPY('All', M, K, A, M, V, M)
         ! Factorize V via dgeqrf
         ! First, query workspace size
         CALL DGEQRF(M, K, V, M, TAU, WORK, NEG_ONE, INFO)
         LWORK = INT(WORK(1))
         DEALLOCATE(WORK)
         ALLOCATE(WORK(LWORK))

         ! Now, factorize V
         CALL DGEQRF(M, K, V, M, TAU, WORK, LWORK, INFO)

         ! Form T. It doesn't matter how we compute it for the first
         ! 3 calls to LARFB{,0C2} But we compute the inverse of our standard
         ! T for the final call to LARFB0C2

         CALL DLARFT('F','C', M, K, V, M, TAU, T, K)
         CALL DLARFT_UT('F','C', '2', M, K, V, M, TAU, T_INV, K)

         CALL CPU_TIME(START)

         CALL DLARFB('L', 'N', 'F', 'C', M, N, K, V, M, T, K, C, M,
     $                  WORK_MAT, M)

         CALL CPU_TIME(FINISH)
         TIMEVEC(1) = FINISH - START
         CALL CPU_TIME(START)

         CALL DLARFB_REF('L', 'N', 'F', 'C', M, N, K, V, M, T, K, C, M,
     $                  WORK_MAT, M)

         CALL CPU_TIME(FINISH)
         TIMEVEC(2) = FINISH - START
         CALL CPU_TIME(START)

         CALL DLARFB0C2('G','B','L','N','F','C',M,N,K,V,M,T,K,C,M)

         CALL CPU_TIME(FINISH)
         TIMEVEC(3) = FINISH - START
         CALL CPU_TIME(START)

         CALL DLARFB0C2('G','A','L','N','F','C',M,N,K,V,M,T,K,C,M)

         CALL CPU_TIME(FINISH)
         TIMEVEC(4) = FINISH - START

         WRITE(*,*) TIMEVEC
10       DEALLOCATE(A)
         DEALLOCATE(V)
         DEALLOCATE(T)
         DEALLOCATE(T_INV)
         DEALLOCATE(C)
         DEALLOCATE(TAU)
         DEALLOCATE(WORK)
      END PROGRAM
