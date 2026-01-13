      PROGRAM TEST_DORGKR
      INTEGER           I, J, K, M, INFO, NEG_ONE, LWORK
      DOUBLE PRECISION  ZERO, ONE, NORM_ORTH, NORM_REPRES, NORM_A
      CHARACTER         JOBT

      DOUBLE PRECISION, ALLOCATABLE :: A(:,:), V(:,:), T(:,:), TAU(:),
     $                  WORK(:), R(:,:), WORKMAT(:,:)
      CHARACTER ::      JOBTS(2)

      EXTERNAL          DORGKR, DLARFT_UT, DLACPY, DLASET

      DOUBLE PRECISION  DLANGE
      EXTERNAL          DLANGE

      JOBTS(1) = '1'
      JOBTS(2) = '2'

      NEG_ONE = -1
      ZERO = 0.0D+0
      ONE = 1.0D+0

      READ(*,*) M, K

      ALLOCATE(A(M,K))
      ALLOCATE(V(M,K))
      ALLOCATE(R(K,K))
      ALLOCATE(T(K,K))
      ALLOCATE(TAU(K))
      ALLOCATE(WORK(1))
      ALLOCATE(WORKMAT(M,M))

      ! Determine necessary workspace sizes
      CALL DGEQRF(M, K, A, M, TAU, WORK, NEG_ONE, INFO)
      LWORK = INT(WORK(1))

      DEALLOCATE(WORK)
      ALLOCATE(WORK(LWORK))

      DO I = 1, 2
         JOBT = JOBTS(I)

         ! Generate A as random data
         CALL RANDOM_NUMBER(A)
         NORM_A = DLANGE('F', M, K, A, M, WORK)
         ! Copy A into V
         CALL DLACPY('All', M, K, A, M, V, M)
         ! Factorize V
         CALL DGEQRF(M, K, V, M, TAU, WORK, LWORK, INFO)
         CALL DLACPY('Upper', K, K, V, M, R, K)
         ! Compute T
         CALL DLARFT('F', 'C', M, K, V, M, TAU, V, M)
         !CALL DLARFT_UT('F', 'C', '2', M, K, V, M, TAU, T, K)
         ! Compute Q
         CALL DORGKR('1', M, K, V, M)

         ! Determine if Q^\top Q - I \approx 0
         CALL DLASET('All', M, M, ZERO, -ONE, WORKMAT, M)
         CALL DSYRK('Upper', 'Transpose', K, M, ONE, V, M, ONE,
     $      WORKMAT, M)

         NORM_ORTH = DLANGE('F', K, K, WORKMAT, M, WORK)
         NORM_ORTH = NORM_ORTH / SQRT(DBLE(K))
         ! Determine if Q*R - A \approx 0
         CALL DLASET('All', M, M, ZERO, ZERO, WORKMAT, M)
         CALL DLACPY('All', M, K, A, M, WORKMAT, M)

         CALL DTRMMOOP('Right', 'Upper', 'N', 'N', 'N', M, K, ONE,
     $      R, K, V, M, -ONE, WORKMAT, M)

         NORM_REPRES = DLANGE('F', M, K, WORKMAT, M, WORK)
         NORM_REPRES = NORM_REPRES / NORM_A

         WRITE(*,*) NORM_ORTH, NORM_REPRES
      END DO


      DEALLOCATE(A)
      DEALLOCATE(V)
      DEALLOCATE(R)
      DEALLOCATE(T)
      DEALLOCATE(TAU)
      DEALLOCATE(WORK)
      DEALLOCATE(WORKMAT)
      END PROGRAM
