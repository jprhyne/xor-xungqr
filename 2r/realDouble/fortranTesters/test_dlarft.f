      SUBROUTINE TEST_DLARFT(M, N)
         ! Arguments
         INTEGER  M, N

         ! Local variables
         DOUBLE PRECISION  NORMA, NORM_ORTH, NORM_REPRES
         INTEGER           LWORK, I, J, INFO
         ! Local arrays
         DOUBLE PRECISION, ALLOCATABLE :: A(:,:), Q(:,:), As(:,:), 
     $            Qs(:,:), WORKMAT(:,:), WORK(:), T(:,:)

         DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TAU

         ! External Subroutines
         EXTERNAL DLACPY, DGEQRF, DLARFT, DTVT, MY_DLARFT

         ! External Functions
         DOUBLE PRECISION, EXTERNAL :: DLANGE

         ! Parameters
         DOUBLE PRECISION ONE, ZERO
         INTEGER           NEG_ONE
         PARAMETER(ONE=1.0D+0, ZERO=0, NEG_ONE=-1.0D+0)

         ALLOCATE(A(M,N))
         ALLOCATE(As(M,N))
         ALLOCATE(Q(M,N))
         ALLOCATE(Qs(M, N))
         ALLOCATE(WORK(1))
         ALLOCATE(TAU(N))
         ALLOCATE(T(N,N))
         ! Generate a random A
         CALL RANDOM_NUMBER(A)

         CALL DLACPY('All', M, N, A, M, As, M)
         NORMA = DLANGE('Frobenius', M, N, A, M, WORK)

         CALL DGEQRF(M, N, A, M, TAU, WORK, NEG_ONE, INFO)
         LWORK = WORK(1)
         DEALLOCATE(WORK)
         ALLOCATE(WORK(LWORK))

         CALL DGEQRF(M, N, A, M, TAU, WORK, LWORK, INFO)
         ! Copy into Q
         CALL DLACPY('All', M, N, A, M, Q, M)
         ! Copy into Qs
         CALL DLACPY('All', M, N, A, M, Qs, M)

         DEALLOCATE(WORK)
         ! Compute the triangular factor T
         !CALL MY_DLARFT_REC(M, N, Q, M, TAU, Q, M)
         CALL MY_DLARFT_UT(M, N, Q, M, TAU, Q, M)

         ! Copy T into where R was inside Q
         !CALL DLACPY('Upper', N, N, T, N, Q, M)

         ! Now call MY_DORG2R
         CALL MY_DORGKR(M, N, Q, M)
         !ALLOCATE(WORK(M*N*M))
         !CALL DORGQR(M, N, N, Qs, M, TAU, WORK, M*N*M, INFO)
         !DEALLOCATE(WORK)
         !ALLOCATE(WORK(N))
         !CALL DORG2R(M,N,N, Qs, M, TAU, WORK, INFO)
         !DEALLOCATE(WORK)
         
         ALLOCATE(WORKMAT(N,N))
         CALL DLASET('A', N,N, ZERO, ZERO, WORKMAT, N)
         CALL DSYRK('Upper', 'Transpose', N, M, ONE, Q, M, ZERO,
     $         WORKMAT, N)

         DO I = 1, N
            WORKMAT(I,I) = WORKMAT(I,I) - 1
         END DO

         NORM_ORTH = DLANGE('Frobenius', N, N, WORKMAT, N, WORK)

         DEALLOCATE(WORKMAT)
         ALLOCATE(WORKMAT(M,N))
         CALL DLACPY('All', M, N, Q, M, WORKMAT, M)
         CALL DTRMM('Right', 'Upper', 'No-transpose', 'non-unit',
     $      M, N, ONE, A, M, WORKMAT, M)

         DO I = 1, M
            DO J = 1, N
               WORKMAT(I,J) = WORKMAT(I,J) - As(I,J)
            END DO
         END DO
         NORM_REPRES = DLANGE('Frobenius', M, N, WORKMAT,
     $      M, WORK)
         NORM_REPRES = NORM_REPRES / NORMA

         DEALLOCATE(WORKMAT)

         WRITE(*,*) "representation norm: ", NORM_REPRES
         WRITE(*,*) "orthogonal norm:     ", NORM_ORTH

         DEALLOCATE(As)
         DEALLOCATE(Q)
         DEALLOCATE(A)
         DEALLOCATE(Qs)
         DEALLOCATE(TAU)
         DEALLOCATE(T)

      END SUBROUTINE
