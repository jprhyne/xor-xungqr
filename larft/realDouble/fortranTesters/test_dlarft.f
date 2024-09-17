      SUBROUTINE TEST_DLARFT(M, N)
         ! Arguments
         INTEGER  M, N

         ! Local variables
         DOUBLE PRECISION  NORMA, NORM_ORTH, NORM_REPRES
         INTEGER           LWORK, I, J, INFO
         CHARACTER         STOREV, DIRECT
         ! Local arrays
         DOUBLE PRECISION, ALLOCATABLE :: A(:,:), Q(:,:), As(:,:), 
     $            Qs(:,:), WORKMAT(:,:), WORK(:), T(:,:), Qt(:,:),
     $            Qts(:,:), Ts(:,:)

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
         ALLOCATE(Qs(M,N))
         ALLOCATE(Qt(N,M))
         ALLOCATE(Qts(N,M))
         ALLOCATE(WORK(1))
         ALLOCATE(TAU(N))
         ALLOCATE(T(N,N))
         ALLOCATE(Ts(N,N))
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
c----------------------------------------------------------------------
         ! Compute the triangular factor T
         DIRECT = 'F'
         STOREV = 'C'
         CALL MY_DLARFT_REC(DIRECT, STOREV, M, N, Q, M, TAU, Q, M)

         ! Test the results
         ! Now call MY_DORGKR
         CALL MY_DORGKR(M, N, Q, M)
         
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

         WRITE(*,*) "Recursive DLARFT. Forward, Column"

         WRITE(*,*) "representation norm: ", NORM_REPRES
         WRITE(*,*) "orthogonal norm:     ", NORM_ORTH

         ! Copy Qs back into Q
         CALL DLACPY('ALL', M, N, Qs, M, Q, M)
c----------------------------------------------------------------------
         ! Above but STOREV = 'R'
         DIRECT = 'F'
         STOREV = 'R'
         ! Copy Qs into Qt
         DO I = 1, M
            DO J = 1, N
               Qt(J,I) = Q(I,J)
            END DO
         END DO
         ! Set T to be all 0s
         DO I = 1, N
            DO J = 1, N
               T(I,J) = ZERO
               Ts(I,J)= ZERO
            END DO
         END DO
         ! Call my dlarft implementation
         CALL MY_DLARFT_REC(DIRECT, STOREV, M, N, Qt, N, TAU, T, N)
         ! Ensure that T is upper triangular
         NORM_REPRES = 0
         DO I = 2, N
            DO J = 1, I-1
               NORM_REPRES = NORM_REPRES + T(I,J) * T(I,J)
               IF(T(I,J).NE.ZERO) THEN
                  WRITE(*,*) "I = ", I, "J = ", J
               END IF
            END DO
         END DO
         IF(NORM_REPRES.NE.ZERO) THEN
            WRITE(*,*) "The upper triangular part of T was touched"
            RETURN
         END IF
         ! Getting to this stage means that T is lower triangular. 
         ! So, we copy T into the upper part of Q
         CALL DLACPY('Upper', N, N, T, N, Q, M)
         ! Next, we check how accurate T was!
         CALL MY_DORGKR(M, N, Q, M)
         
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

         WRITE(*,*) "Recursive DLARFT. Forward, Row"

         WRITE(*,*) "representation norm: ", NORM_REPRES
         WRITE(*,*) "orthogonal norm:     ", NORM_ORTH

         ! Copy Qs back into Q
         CALL DLACPY('ALL', M, N, Qs, M, Q, M)

c----------------------------------------------------------------------
         ! Copy Qs back into Q
         CALL DLACPY('ALL', M, N, Qs, M, Q, M)

         ! Compute the triangular factor T
         CALL MY_DLARFT_UT(M, N, Q, M, TAU, Q, M)

         ! Now call MY_DORGKR
         CALL MY_DORGKR(M, N, Q, M)
         
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

         WRITE(*,*) "UT DLARFT"

         WRITE(*,*) "representation norm: ", NORM_REPRES
         WRITE(*,*) "orthogonal norm:     ", NORM_ORTH

         DEALLOCATE(As)
         DEALLOCATE(Q)
         DEALLOCATE(A)
         DEALLOCATE(Qs)
         DEALLOCATE(TAU)
         DEALLOCATE(T)

      END SUBROUTINE
