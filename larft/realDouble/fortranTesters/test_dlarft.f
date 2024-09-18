      SUBROUTINE TEST_DLARFT(M, N)
         ! Arguments
         INTEGER  M, N

         ! Local variables
         DOUBLE PRECISION  NORMA, NORM_ORTH, NORM_REPRES, TMP
         INTEGER           LWORK, I, J, INFO
         CHARACTER         STOREV, DIRECT
         ! Local arrays
         DOUBLE PRECISION, ALLOCATABLE :: A(:,:), Q(:,:), As(:,:), 
     $            Qs(:,:), WORKMAT(:,:), WORK(:), T(:,:), Qt(:,:),
     $            Qts(:,:), Ts(:,:), At(:,:), Ats(:,:)

         DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TAU

         ! External Subroutines
         EXTERNAL DLACPY, DGEQRF, DLARFT, DTVT, MY_DLARFT

         ! External Functions
         DOUBLE PRECISION, EXTERNAL :: DLANGE

         ! Parameters
         DOUBLE PRECISION ONE, ZERO
         INTEGER           NEG_ONE
         PARAMETER(ONE=1.0D+0, ZERO=0, NEG_ONE=-1)

         ALLOCATE(A(M,N))
         ALLOCATE(At(N,M))
         ALLOCATE(Ats(N,M))
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
c----------------------------------------------------------------------
         ! A = QR
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
c----------------------------------------------------------------------
         ! Above but STOREV = 'R'
         DIRECT = 'F'
         STOREV = 'R'
         ! A = LQ
         ! Copy As into At
         DO I = 1, M
            DO J = 1, N
               At(J,I) = As(I,J)
            END DO
         END DO
         ALLOCATE(WORK(1))
         ! Determine how much workspace is needed
         CALL DGELQF(N, M, At, N, TAU, WORK, NEG_ONE, INFO)
         LWORK = WORK(1)
         DEALLOCATE(WORK)
         ALLOCATE(WORK(LWORK))
         CALL DGELQF(N, M, At, N, TAU, WORK, LWORK, INFO)
         DEALLOCATE(WORK)
         ! Now, we have the reflectors for Q. Store this in Ats in order
         ! to store L and ensure we don't modify the reflectors on exit
         CALL DLACPY('All', N, M, At, N, Ats, N)
         ! Set T to be all 0s
         DO I = 1, N
            DO J = 1, N
               T(I,J) = ZERO
            END DO
         END DO
         ! Call my dlarft implementation
         CALL MY_DLARFT_REC(DIRECT, STOREV, M, N, At, N, TAU, T, N)
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
            WRITE(*,*) "The lower triangular part of T was touched"
            RETURN
         END IF
         ! Make sure that At was not touched
         NORM_REPRES = 0.0D+0
         DO J = 1, M
            DO I = 1, N
               TMP = Ats(I,J) - At(I,J)
               NORM_REPRES = NORM_REPRES + TMP*TMP
               IF(TMP.NE.ZERO) THEN
                  WRITE(*,*) "I = ", I, "J = ", J
               END IF
            END DO
         END DO
         IF(NORM_REPRES.NE.ZERO) THEN
            WRITE(*,*) "A was changed"
            RETURN
         END IF
         ! Getting to this stage means that T is upper triangular. 
         ! We must now construct Q in order to check if we have a correct T
         ! Right now, I only have an implementation for my_dorgkr for when
         ! we want to compute the Q from the QR factorization. So, we will 
         ! Copy over the transpose of the reflectors into Q
         DO I = 1, N
            DO J = I+1, M
               Q(J,I) = At(I,J)
            END DO
         END DO
         ! Copy in the T matrix to the upper triangular part of Q
         DO I = 1, N
            DO J = I, N
               Q(I,J) = T(I,J)
            END DO
         END DO
         ! Now we can call DORGKR
         CALL MY_DORGKR(M, N, Q, M)
         ! Next, we test that Q is orthogonal
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
         ! Now, we test that A = Q(L^\top) = QR
         CALL DLACPY('All', M, N, Q, M, WORKMAT, M)
         CALL DTRMM('Right', 'Lower', 'Transpose', 'Non-unit',
     $      M, N, ONE, At, N, WORKMAT, M)

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
         RETURN
c----------------------------------------------------------------------
         ! Above but DIRECT = 'B' and STOREV = 'C'
         ! A = QL
         DIRECT = 'B'
         STOREV = 'C'
         ! Copy As into A
         CALL DLACPY('All', M, N, As, M, A, M)

         ALLOCATE(WORK(1))
         CALL DGEQLF(M, N, A, M, TAU, WORK, NEG_ONE, INFO)
         LWORK = WORK(1)
         DEALLOCATE(WORK)
         ALLOCATE(WORK(LWORK))
         CALL DGEQLF(M, N, A, M, TAU, WORK, LWORK, INFO)
         ! Set T to be all 0s
         DO I = 1, N
            DO J = 1, N
               T(I,J) = ZERO
               Ts(I,J) = ZERO
            END DO
         END DO
         ! Call my dlarft implementation
         CALL MY_DLARFT_REC(DIRECT, STOREV, M, N, A, M, TAU, T, N)
         ! Ensure that A was not touched
         ! Since we don't have a nice helper (yet?) to compute the Q associated
         ! with the QL that takes in a T matrix, we test the accuracy of T by
         ! comparing it against the result that DLARFT_REF will give
         !
         ! Call reference dlarft
         CALL DLARFT(DIRECT, STOREV, M, N, A, M, TAU, Ts, N)
         
         ! Compute the forward error
         NORM_REPRES = 0.0D+0
         DO I = 1, N
            DO J = 1, N
               TMP = Ts(I,J) - T(I,J)
               NORM_REPRES = NORM_REPRES + TMP*TMP
            END DO
         END DO
         ! Compute NormT
         TMP = 0.0D+0
         DO I = 1, N
            DO J = 1, N
               TMP = TMP + Ts(I,J)*Ts(I,J)
            END DO
         END DO
         NORM_REPRES = SQRT(NORM_REPRES)/TMP

         WRITE(*,*) "Recursive DLARFT. Backwards, Col"

         WRITE(*,*) "|T_{mine} - T_{ref}|_F/|T_{ref}|_F: ", NORM_REPRES
c----------------------------------------------------------------------
         ! Above but STOREV = 'R'
         ! A = RQ
         DIRECT = 'B'
         STOREV = 'R'
         ! Copy As into At
         DO I = 1, M
            DO J = 1, N
               At(J,I) = As(I,J)
            END DO
         END DO
         ALLOCATE(WORK(1))
         CALL DGERQF(N, M, At, N, TAU, WORK, NEG_ONE, INFO)
         LWORK = WORK(1)
         DEALLOCATE(WORK)
         ALLOCATE(WORK(LWORK))
         CALL DGERQF(N, M, At, N, TAU, WORK, LWORK, INFO)
         ! Set T to be all 0s
         DO I = 1, N
            DO J = 1, N
               T(I,J) = ZERO
               Ts(I,J)= ZERO
            END DO
         END DO
         ! Call my dlarft implementation
         CALL MY_DLARFT_REC(DIRECT, STOREV, M, N, A, M, TAU, T, N)
         ! Ensure that A was not touched
         ! Since we don't have a nice helper (yet?) to compute the Q associated
         ! with the QL that takes in a T matrix, we test the accuracy of T by
         ! comparing it against the result that DLARFT_REF will give
         !
         ! Call reference dlarft
         CALL DLARFT(DIRECT, STOREV, M, N, A, M, TAU, Ts, N)
         
         ! Compute the forward error
         NORM_REPRES = 0.0D+0
         DO I = 1, N
            DO J = 1, N
               TMP = Ts(I,J) - T(I,J)
               NORM_REPRES = NORM_REPRES + TMP*TMP
            END DO
         END DO
         ! Compute NormT
         TMP = 0.0D+0
         DO I = 1, N
            DO J = 1, N
               TMP = TMP + Ts(I,J)*Ts(I,J)
            END DO
         END DO
         NORM_REPRES = SQRT(NORM_REPRES)/TMP

         WRITE(*,*) "Recursive DLARFT. Backwards, Row"

         WRITE(*,*) "|T_{mine} - T_{ref}|_F/|T_{ref}|_F: ", NORM_REPRES

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
         DEALLOCATE(At)
         DEALLOCATE(Ats)
         DEALLOCATE(At)
         DEALLOCATE(Qs)
         DEALLOCATE(TAU)
         DEALLOCATE(T)
         DEALLOCATE(Ts)

      END SUBROUTINE
