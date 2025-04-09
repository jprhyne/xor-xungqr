      SUBROUTINE TEST_DGEQRF(M, N)
         ! Arguments
         INTEGER  M, N

         ! Local variables
         DOUBLE PRECISION  NORM_FORWARD, TMP, EPS, NORM_A, DN
         INTEGER           I, LWORK, INFO
         CHARACTER         STOREV, DIRECT
         ! Local arrays
         DOUBLE PRECISION, ALLOCATABLE :: A(:,:), As(:,:),
     $            WORKMAT(:,:), WORK(:), TAU(:), R(:,:)

         ! Intrinsic subroutines
         INTRINSIC SQRT 
         ! External Subroutines
         EXTERNAL DLACPY, DGEQRF_QR3, DORGQR

         ! External Functions
         DOUBLE PRECISION, EXTERNAL :: DLANGE

         ! Parameters
         DOUBLE PRECISION ONE, ZERO, MONE
         INTEGER           NEG_ONE
         PARAMETER(ONE=1.0D+0, ZERO=0.0D+0, MONE=-1.0D+0, NEG_ONE=-1)

         DN = N
         
         ALLOCATE(A(M,N))
         ALLOCATE(As(M,N))
         ALLOCATE(R(N,N))
         ALLOCATE(TAU(N))
         ALLOCATE(WORKMAT(M,M))
         EPS = EPSILON(ONE)
         ! Generate a random A
         CALL RANDOM_NUMBER(A)
         ! Compute a norm of A
         NORM_A = DLANGE('Frobenius', M, N, A, M, WORK)
         ! Store a copy of A
         CALL DLACPY('All', M, N, A, M, As, M)
c----------------------------------------------------------------------
         ! A = QR
         ! Workspace query
         ALLOCATE(WORK(1))
         CALL DGEQRF_QR3(M, N, 32, A, M, TAU, WORK, NEG_ONE, INFO)
         LWORK = WORK(1)
         ! Query DORGQR
         CALL DORGQR(M, N, N, A, M, TAU, WORK, NEG_ONE, INFO)
         IF (LWORK < WORK(1)) THEN
            LWORK = WORK(1)
         END IF
         DEALLOCATE(WORK)
         ALLOCATE(WORK(LWORK))

         ! Set R to be 0
         CALL DLASET('All', N, N, ZERO, ZERO, R, N)

         ! Compute the Q factor
         CALL DGEQRF_QR3(M, N, 32, A, M, TAU, WORK, LWORK, INFO)
         ! Copy in R
         CALL DLACPY('Upper', N, N, A, M, R, N)
         ! Compute Q
         CALL DORGQR(M, N, N, A, M, TAU, WORK, LWORK, INFO)

         ! Compute ||QR - A||
         ! Copy A into WORKMAT
         CALL DLACPY('All', M, N, As, M, WORKMAT, M)
         ! Compute Workmat = QR - A
         CALL DGEMM('No transpose', 'No transpose', M, N, N, ONE,
     $               A, M, R, N, MONE, WORKMAT, M)
         ! Compute the norm
         NORM_FORWARD = DLANGE('Frobenius', M, N, WORKMAT, M, WORK)
         ! Divide by the norm of A
         NORM_FORWARD = NORM_FORWARD / NORM_A


         WRITE(*,*) "Representation Error", NORM_FORWARD
         ! WORKMAT = Q**T Q
         CALL DGEMM('Transpose', 'No transpose', N, N, M, ONE, A, M,
     $               A, M, ZERO, WORKMAT, M)
         ! Compute WORKMAT = WORKMAT - I
         DO I = 1, N
            WORKMAT(I,I) = WORKMAT(I,I) - ONE
         END DO
         ! Compute the norm
         NORM_FORWARD = DLANGE('Frobenius', N, N, WORKMAT, M, WORK)
         ! Divide by the norm of I
         NORM_FORWARD = NORM_FORWARD / SQRT(DN)
         WRITE(*,*) "Orthogonal Error", NORM_FORWARD


         DEALLOCATE(WORK)
   10    DEALLOCATE(A)
         DEALLOCATE(As)
         DEALLOCATE(TAU)

      END SUBROUTINE
