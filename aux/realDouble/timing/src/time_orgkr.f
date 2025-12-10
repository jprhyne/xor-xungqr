      PROGRAM TIME_ORGKR
         INTEGER M,N,LWORK,INFO,NEG_ONE
         REAL START, FINISH

         DOUBLE PRECISION, ALLOCATABLE :: A(:,:), WORK(:), TAU(:),
     $      V(:,:), T(:,:)

         DOUBLE PRECISION, DIMENSION(3) :: TIMES

         INTRINSIC CPU_TIME, MAX

         READ(*,*) M,N

         LWORK = MAX(M,N)
         NEG_ONE = -1

         ALLOCATE(A(M,N))
         ALLOCATE(T(N,N))
         ALLOCATE(V(M,N))
         ALLOCATE(WORK(LWORK))
         ALLOCATE(TAU(LWORK))

         CALL DGEQRF(M, N, A, M, TAU, WORK, NEG_ONE, INFO)
         LWORK = MAX(INT(WORK(1)), LWORK)

         DEALLOCATE(WORK)
         ALLOCATE(WORK(LWORK))

         CALL RANDOM_NUMBER(A)

         CALL DGEQRF(M, N, A, M, TAU, WORK, LWORK, INFO)
         CALL DLACPY('All', M, N, A, M, V, M)

         CALL CPU_TIME(START)
         CALL DORG2R(M, N, N, V, M, TAU, WORK, INFO)
         CALL CPU_TIME(FINISH)
         TIMES(1) = FINISH - START

         CALL DLACPY('All', M, N, A, M, V, M)

         CALL CPU_TIME(START)
         CALL DORG2R_REF(M, N, N, V, M, TAU, WORK, INFO)
         CALL CPU_TIME(FINISH)
         TIMES(2) = FINISH - START

         CALL DLACPY('All', M, N, A, M, V, M)

         CALL DLARFT('F', 'C', M, N, V, M, TAU, T, N)
         CALL DLACPY('Upper', N, N, T, N, V, M)

         CALL CPU_TIME(START)
         CALL DORGKR(M, N, A, M)
         CALL CPU_TIME(FINISH)
         TIMES(3) = FINISH - START

         WRITE(*,*) "QR_TIMES: ", TIMES

         DEALLOCATE(A)
         DEALLOCATE(T)
         DEALLOCATE(V)
         DEALLOCATE(WORK)
         DEALLOCATE(TAU)
         !
         ALLOCATE(A(N,M))
         ALLOCATE(T(N,N))
         ALLOCATE(V(N,M))
         ALLOCATE(WORK(LWORK))
         ALLOCATE(TAU(N))

         CALL DGELQF(N, M, A, N, TAU, WORK, NEG_ONE, INFO)
         LWORK = MAX(INT(WORK(1)), LWORK)

         DEALLOCATE(WORK)
         ALLOCATE(WORK(LWORK))

         CALL RANDOM_NUMBER(A)

         CALL DGELQF(N, M, A, N, TAU, WORK, LWORK, INFO)
         CALL DLACPY('All', N, M, A, N, V, N)

         CALL CPU_TIME(START)
         CALL DORGL2(N, M, N, V, N, TAU, WORK, INFO)
         CALL CPU_TIME(FINISH)
         TIMES(1) = FINISH - START

         CALL DLACPY('All', M, N, A, M, V, M)

         CALL CPU_TIME(START)
         CALL DORGL2_REF(N, M, N, V, N, TAU, WORK, INFO)
         CALL CPU_TIME(FINISH)
         TIMES(2) = FINISH - START

         CALL DLACPY('All', M, N, A, M, V, M)

         CALL DLARFT('F', 'T', M, N, V, N, TAU, T, N)
         CALL DLACPY('Lower', N, N, T, N, V, N)

         CALL CPU_TIME(START)
         CALL DORGLK(N, M, A, N)
         CALL CPU_TIME(FINISH)
         TIMES(3) = FINISH - START

         WRITE(*,*) "LQ_TIMES: ", TIMES
         DEALLOCATE(A)
         DEALLOCATE(T)
         DEALLOCATE(V)
         DEALLOCATE(WORK)
         DEALLOCATE(TAU)
      END PROGRAM
