      SUBROUTINE TEST_DORGKR(M, N, K, LDA, LDQ)
         ! Arguments
         INTEGER  M, N, K, LDA, LDQ

         ! Local variables
         DOUBLE PRECISION  NORMA, NORM_ORTH, NORM_REPRES
         INTEGER           LWORK, I, J, INFO
         ! Local arrays
         DOUBLE PRECISION, ALLOCATABLE :: A(:,:), Q(:,:), As(:,:), 
     $            Qs(:,:), WORKMAT(:,:), WORK(:), T(:,:)

         DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TAU

         ! External Subroutines
         EXTERNAL DLACPY, DGEQRF, DLARFT, DTVT

         ! External Functions
         DOUBLE PRECISION, EXTERNAL :: DLANGE

         ! Parameters
         DOUBLE PRECISION ONE, ZERO
         INTEGER           NEG_ONE
         PARAMETER(ONE=1.0D+0, ZERO=0, NEG_ONE=-1.0D+0)

         IF(LDQ.LT.M.OR.LDA.LT.M) THEN
            WRITE(*,*) "Incorrect LDA or LDQ value"
            RETURN
         ELSE IF(K.GT.N) THEN
            WRITE(*,*) "Incorrect K value"
         END IF

         ALLOCATE(A(LDA,K))
         ALLOCATE(As(LDA,K))
         ALLOCATE(Q(LDQ,K))
         ALLOCATE(Qs(LDQ, N))
         ALLOCATE(WORK(1))
         ALLOCATE(TAU(K))
         ALLOCATE(T(K,K))
         ! Generate a random A
         CALL RANDOM_NUMBER(A)

         CALL DLACPY('All', M, K, A, LDA, As, LDA)
         NORMA = DLANGE('Frobenius', M, K, A, LDA, WORK)

         CALL DGEQRF(M, K, A, LDA, TAU, WORK, NEG_ONE, INFO)
         LWORK = WORK(1)
         DEALLOCATE(WORK)
         ALLOCATE(WORK(LWORK))

         CALL DGEQRF(M, K, A, LDA, TAU, WORK, LWORK, INFO)
         ! Copy into Q
         CALL DLACPY('All', M, K, A, LDA, Q, LDQ)
         ! Copy into Qs
         CALL DLACPY('All', M, K, A, LDA, Qs, LDQ)

         DEALLOCATE(WORK)
         ! Compute the triangular factor T
         !CALL DLARFT('Forward', 'Column', M, K, Q, LDQ, TAU, Q, LDQ)
         CALL MY_DLARFT_REC(M, K, Q, LDQ, TAU, Q, LDQ)
         !CALL MY_DLARFT_UT(M, K, Q, LDQ, TAU, Q, LDQ)

         ! Copy T into where R was inside Q
         !CALL DLACPY('Upper', K, K, T, K, Q, LDQ)

         ! Now call MY_DORG2R
         CALL MY_DORGKR(M, K, Q, LDQ)
         !ALLOCATE(WORK(M*N*M))
         !CALL DORGQR(M, N, K, Qs, LDQ, TAU, WORK, M*N*M, INFO)
         !DEALLOCATE(WORK)
         !ALLOCATE(WORK(N))
         !CALL DORG2R(M,N,K, Qs, LDQ, TAU, WORK, INFO)
         !DEALLOCATE(WORK)
         
         ALLOCATE(WORKMAT(K,K))
         CALL DLASET('A', K,K, ZERO, ZERO, WORKMAT, K)
         CALL DSYRK('Upper', 'Transpose', K, M, ONE, Q, LDQ, ZERO,
     $         WORKMAT, K)

         DO I = 1, K
            WORKMAT(I,I) = WORKMAT(I,I) - 1
         END DO

         NORM_ORTH = DLANGE('Frobenius', K, K, WORKMAT, K, WORK)

         DEALLOCATE(WORKMAT)
         ALLOCATE(WORKMAT(M,K))
         CALL DLACPY('All', M, K, Q, LDQ, WORKMAT, M)
         CALL DTRMM('Right', 'Upper', 'No-transpose', 'non-unit',
     $      M, K, ONE, A, LDA, WORKMAT, M)

         DO I = 1, M
            DO J = 1, K
               WORKMAT(I,J) = WORKMAT(I,J) - As(I,J)
            END DO
         END DO
         NORM_REPRES = DLANGE('Frobenius', M, K, WORKMAT,
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
