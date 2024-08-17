      SUBROUTINE TEST_ZLARFT(M, N)
         ! Arguments
         INTEGER  M, N

         ! Local variables
         DOUBLE PRECISION  NORMA, NORM_ORTH, NORM_REPRES
         DOUBLE PRECISION  TMP1, TMP2
         INTEGER           LWORK, I, J, INFO
         ! Local arrays
         COMPLEX*16, ALLOCATABLE :: A(:,:), Q(:,:), As(:,:), 
     $            Qs(:,:), WORKMAT(:,:), WORK(:), T(:,:)

         COMPLEX*16, DIMENSION(:), ALLOCATABLE :: TAU

         ! External Subroutines
         EXTERNAL ZLACPY, ZGEQRF, ZLARFT, ZTVT, MY_ZLARFT

         ! External Functions
         DOUBLE PRECISION, EXTERNAL :: ZLANGE

         ! Parameters
         DOUBLE PRECISION ONE, ZERO
         COMPLEX*16       CONE,CZERO
         INTEGER          NEG_ONE
         PARAMETER(ONE=1.0D+0, ZERO=0, NEG_ONE=-1.0D+0, 
     $            CONE=COMPLEX(1,0), CZERO=COMPLEX(0,0))

         ALLOCATE(A(M,N))
         ALLOCATE(As(M,N))
         ALLOCATE(Q(M,N))
         ALLOCATE(Qs(M, N))
         ALLOCATE(WORK(1))
         ALLOCATE(TAU(N))
         ALLOCATE(T(N,N))
         ! Generate a random complex A
         DO I = 1, M
            DO J = 1, N
               CALL RANDOM_NUMBER(TMP1)
               CALL RANDOM_NUMBER(TMP2)
               A(I,J) = COMPLEX(TMP1, TMP2)
            END DO
         END DO

         ! Copy A into As
         CALL ZLACPY('All', M, N, A, M, As, M)

         ! Store ||A||_F for later use
         NORMA = ZLANGE('Frobenius', M, N, A, M, WORK)

         ! Compute the QR decomposition of A

         ! Workspace query
         CALL ZGEQRF(M, N, A, M, TAU, WORK, NEG_ONE, INFO)
         LWORK = WORK(1)
         DEALLOCATE(WORK)
         ALLOCATE(WORK(LWORK))

         ! Doing the factorization
         CALL ZGEQRF(M, N, A, M, TAU, WORK, LWORK, INFO)
         ! Copy into Q
         CALL ZLACPY('All', M, N, A, M, Q, M)
         ! Copy into Qs
         CALL ZLACPY('All', M, N, A, M, Qs, M)

         DEALLOCATE(WORK)
         ! Compute the triangular factor T
         CALL MY_ZLARFT_UT_V2(M, N, Q, M, TAU, Q, M)

         ! Now call MY_ZORGKR
         CALL MY_ZORGKR(M, N, Q, M)
         
         ALLOCATE(WORKMAT(N,N))
         ! Compute W = Q**T *Q - I
         ! W = 0
         CALL ZLASET('A', N,N, ZERO, ZERO, WORKMAT, N)
         ! W = Q**T * Q
         CALL ZSYRK('Upper', 'Transpose', N, M, ONE, Q, M, ZERO,
     $         WORKMAT, N)
         ! Q = W - I = Q**T * Q - I
         DO I = 1, N
            WORKMAT(I,I) = WORKMAT(I,I) - 1
         END DO

         ! Compute ||Q**T * Q - I||_F
         NORM_ORTH = ZLANGE('Frobenius', N, N, WORKMAT, N, WORK)

         DEALLOCATE(WORKMAT)
         ALLOCATE(WORKMAT(M,N))
         ! Compute ||A - Q*R||_F / ||A||_F
         ! W = Q
         CALL ZLACPY('All', M, N, Q, M, WORKMAT, M)
         ! W = W*R = Q*R
         CALL ZTRMM('Right', 'Upper', 'No-transpose', 'Non-unit',
     $      M, N, ONE, A, M, WORKMAT, M)

         ! W = W - A = QR - A
         DO I = 1, M
            DO J = 1, N
               WORKMAT(I,J) = WORKMAT(I,J) - As(I,J)
            END DO
         END DO
         ! NORM_REPRES = ||A - Q*R||_F / ||A||_F
         NORM_REPRES = ZLANGE('Frobenius', M, N, WORKMAT,
     $      M, WORK)
         NORM_REPRES = NORM_REPRES / NORMA

         DEALLOCATE(WORKMAT)

         ! Print out the information to the console
         WRITE(*,*) "representation norm: ", NORM_REPRES
         WRITE(*,*) "orthogonal norm:     ", NORM_ORTH

         ! Free our memory
         DEALLOCATE(As)
         DEALLOCATE(Q)
         DEALLOCATE(A)
         DEALLOCATE(Qs)
         DEALLOCATE(TAU)
         DEALLOCATE(T)

      END SUBROUTINE
