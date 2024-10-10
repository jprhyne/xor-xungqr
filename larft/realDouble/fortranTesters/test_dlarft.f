      SUBROUTINE TEST_DLARFT(M, N)
         ! Arguments
         INTEGER  M, N

         ! Local variables
         DOUBLE PRECISION  NORMT, NORM_FORWARD, TMP, EPS
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
         PARAMETER(ONE=1.0D+0, ZERO=0.0D+0, NEG_ONE=-1)

         ALLOCATE(A(M,N))
         ALLOCATE(At(N,M))
         ALLOCATE(Ats(N,M))
         ALLOCATE(As(M,N))
         ALLOCATE(Q(M,N))
         ALLOCATE(Qs(M,N))
         ALLOCATE(Qt(N,M))
         ALLOCATE(Qts(N,M))
         ALLOCATE(TAU(N))
         ALLOCATE(T(N,N))
         ALLOCATE(Ts(N,N))
         EPS = EPSILON(ONE)
         ! Generate a random A
         CALL RANDOM_NUMBER(A)

         CALL DLACPY('All', M, N, A, M, As, M)
c----------------------------------------------------------------------
         ! A = QR
         ALLOCATE(WORK(1))
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
         ! Set T and Ts to be 0
         DO I = 1, N
            DO J = 1, N
               T(I,J)  = ZERO
               Ts(I,J) = ZERO
            END DO
         END DO
         ! Compute the triangular factor T
         DIRECT = 'F'
         STOREV = 'C'
         CALL MY_DLARFT_REC(DIRECT, STOREV, M, N, Q, M, TAU, T, N)
         ! Ensure that T is upper triangular
         NORM_FORWARD = 0.0D+0
         DO I = 2, N
            DO J = 1, I-1
               NORM_FORWARD = NORM_FORWARD + T(I,J) * T(I,J)
               IF(T(I,J).NE.ZERO) THEN
                  WRITE(*,*) "I = ", I, "J = ", J
               END IF
            END DO
         END DO
         IF(NORM_FORWARD.NE.ZERO) THEN
            WRITE(*,*) "The lower triangular part of T was touched"
            GOTO 10
         END IF
         ! Make sure that Q was not touched
         NORM_FORWARD = 0.0D+0
         DO I = 1, M
            DO J = 1, N
               TMP = Q(I,J) - Qs(I,J)
               NORM_FORWARD = NORM_FORWARD + TMP*TMP
               IF(TMP.NE.ZERO) THEN
                  WRITE(*,*) "I = ", I, "J = ", J
               END IF
            END DO
         END DO
         IF(NORM_FORWARD.NE.ZERO) THEN
            WRITE(*,*) "A was changed"
            GOTO 10
         END IF

         ! We compare the result to the existing dlarft reference implementation
         CALL DLARFT_REF(DIRECT, STOREV, M, N, Q, M, TAU, Ts, N)
         ! Test the results
         NORM_FORWARD = 0.0D+0
         NORMT = DLANGE('Frobenius', N, N, Ts, N, WORK)

         DO I = 1, N
            DO J = 1, N
               TMP = Ts(I,J) - T(I,J)
               TMP = TMP * TMP
               NORM_FORWARD = NORM_FORWARD + TMP
               IF(TMP.GT.EPS*NORMT) THEN
                  WRITE(*,*) "I: ", I, "J: ", J, "Val: ", TMP
               END IF
            END DO
         END DO

         NORM_FORWARD = SQRT(NORM_FORWARD)

         WRITE(*,*) "Recursive DLARFT. Forward, Column"

         WRITE(*,*) "Forward Error with Reference DLARFT: ", 
     $               NORM_FORWARD/NORMT
c----------------------------------------------------------------------
         ! Above but STOREV = 'R'
         ! Set T and Ts to be 0
         DO I = 1, N
            DO J = 1, N
               T(I,J)  = ZERO
               Ts(I,J) = ZERO
            END DO
         END DO
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
         ! Call my dlarft implementation
         CALL MY_DLARFT_REC(DIRECT, STOREV, M, N, At, N, TAU, T, N)
         ! Ensure that T is upper triangular
         NORM_FORWARD = 0
         DO I = 2, N
            DO J = 1, I-1
               NORM_FORWARD = NORM_FORWARD + T(I,J) * T(I,J)
               IF(T(I,J).NE.ZERO) THEN
                  WRITE(*,*) "I = ", I, "J = ", J
               END IF
            END DO
         END DO
         IF(NORM_FORWARD.NE.ZERO) THEN
            WRITE(*,*) "The lower triangular part of T was touched"
            GOTO 10
         END IF
         ! Make sure that At was not touched
         NORM_FORWARD = 0.0D+0
         DO J = 1, M
            DO I = 1, N
               TMP = Ats(I,J) - At(I,J)
               NORM_FORWARD = NORM_FORWARD + TMP*TMP
               IF(TMP.NE.ZERO) THEN
                  WRITE(*,*) "I = ", I, "J = ", J
               END IF
            END DO
         END DO
         IF(NORM_FORWARD.NE.ZERO) THEN
            WRITE(*,*) "A was changed"
            GOTO 10
         END IF
         ! Getting to this stage means that T is upper triangular. 
         ! Now, we compare our result against dlarft_ref
         CALL DLARFT_REF(DIRECT, STOREV, M, N, At, N, TAU, Ts, N)
         ! Test the results
         NORM_FORWARD = 0.0D+0
         NORMT = DLANGE('Frobenius', N, N, Ts, N, WORK)

         DO I = 1, N
            DO J = 1, N
               TMP = Ts(I,J) - T(I,J)
               TMP = TMP * TMP
               NORM_FORWARD = NORM_FORWARD + TMP
               IF(TMP.GT.EPS*NORMT) THEN
                  WRITE(*,*) "I: ", I, "J: ", J, "Val: ", TMP
               END IF
            END DO
         END DO

         NORM_FORWARD = SQRT(NORM_FORWARD)
         WRITE(*,*) "Recursive DLARFT. Forward, Row"

         WRITE(*,*) "Forward Error with Reference DLARFT: ", 
     $               NORM_FORWARD/NORMT
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
         DEALLOCATE(WORK)
         ! Copy into Q
         CALL DLACPY('All', M, N, A, M, Q, M)
         ! Set T to be all 0s
         DO I = 1, N
            DO J = 1, N
               T(I,J) = ZERO
               Ts(I,J) = ZERO
            END DO
         END DO
         ! Call my dlarft implementation
         CALL MY_DLARFT_REC(DIRECT, STOREV, M, N, Q, M, TAU, T, N)
         ! Ensure that T is lower triangular
         NORM_FORWARD = 0.0D+0
         DO I = 1, N-1
            DO J = I+1, N
               NORM_FORWARD = NORM_FORWARD + T(I,J) * T(I,J)
               IF(T(I,J).NE.ZERO) THEN
                  WRITE(*,*) "I = ", I, "J = ", J
               END IF
            END DO
         END DO
         IF(NORM_FORWARD.NE.ZERO) THEN
            WRITE(*,*) "The upper triangular part of T was touched"
            GOTO 10
         END IF
         ! Make sure that Q was not touched
         NORM_FORWARD = 0.0D+0
         DO I = 1, M
            DO J = 1, N
               TMP = Q(I,J) - A(I,J)
               NORM_FORWARD = NORM_FORWARD + TMP*TMP
               IF(TMP.NE.ZERO) THEN
                  WRITE(*,*) "I = ", I, "J = ", J
               END IF
            END DO
         END DO
         IF(NORM_FORWARD.NE.ZERO) THEN
            WRITE(*,*) "A was changed"
            GOTO 10
         END IF
         ! Now, we compare our result against dlarft_rec
         CALL DLARFT_REF(DIRECT, STOREV, M, N, Q, M, TAU, Ts, N)
         ! Test the results
         NORM_FORWARD = 0.0D+0
         NORMT = DLANGE('Frobenius', N, N, Ts, N, WORK)

         DO I = 1, N
            DO J = 1, N
               TMP = Ts(I,J) - T(I,J)
               NORM_FORWARD = NORM_FORWARD + TMP * TMP
               IF(TMP.GT.EPS*NORMT) THEN
                  WRITE(*,*) "I: ", I, "J: ", J, "Val: ", TMP
               END IF
            END DO
         END DO

         NORM_FORWARD = SQRT(NORM_FORWARD)

         WRITE(*,*) "Recursive DLARFT. Backwards, Col"

         WRITE(*,*) "Forward Error with Reference DLARFT: ", 
     $               NORM_FORWARD/NORMT
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
         DEALLOCATE(WORK)
         CALL DLACPY('All', N, M, At, N, Qt, N)
         ! Set T to be all 0s
         DO I = 1, N
            DO J = 1, N
               T(I,J) = ZERO
               Ts(I,J)= ZERO
            END DO
         END DO
         ! Call my dlarft implementation
         CALL MY_DLARFT_REC(DIRECT, STOREV, M, N, Qt, N, TAU, T, N)
         ! Ensure that T is lower triangular
         NORM_FORWARD = 0.0D+0
         DO I = 1, N-1
            DO J = I+1, N
               NORM_FORWARD = NORM_FORWARD + T(I,J) * T(I,J)
               IF(T(I,J).NE.ZERO) THEN
                  WRITE(*,*) "I = ", I, "J = ", J
               END IF
            END DO
         END DO
         IF(NORM_FORWARD.NE.ZERO) THEN
            WRITE(*,*) "The upper triangular part of T was touched"
            GOTO 10
         END IF
         ! Ensure that Qt was not touched
         DO J = 1, M
            DO I = 1, N
               TMP = Qt(I, J) - At(I, J)
               TMP = TMP*TMP
               IF(TMP.NE.0) THEN
                  WRITE(*,*) "A was modified"
                  GOTO 10
               END IF
            END DO
         END DO
         ! Call reference dlarft
         CALL DLARFT_REF(DIRECT, STOREV, M, N, At, N, TAU, Ts, N)
         
         NORMT = DLANGE('Frobenius', N, N, Ts, N, WORK)
         ! Compute the forward error
         NORM_FORWARD = 0.0D+0
         DO I = 1, N
            DO J = 1, N
               TMP = Ts(I,J) - T(I,J)
               NORM_FORWARD = NORM_FORWARD + TMP*TMP
               IF(TMP.GT.NORMT*EPS) THEN
                  WRITE(*,*) "I = ", I, "J = ", J
               END IF
            END DO
         END DO
         NORM_FORWARD = SQRT(NORM_FORWARD)/NORMT

         WRITE(*,*) "Recursive DLARFT. Backwards, Row"

         WRITE(*,*) "|T_{mine} - T_{ref}|_F/|T_{ref}|_F: ", NORM_FORWARD
         GOTO 10

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

         DEALLOCATE(WORKMAT)


         WRITE(*,*) "representation norm: ", NORM_FORWARD

   10    DEALLOCATE(A)
         DEALLOCATE(At)
         DEALLOCATE(Ats)
         DEALLOCATE(As)
         DEALLOCATE(Q)
         DEALLOCATE(Qs)
         DEALLOCATE(Qt)
         DEALLOCATE(Qts)
         DEALLOCATE(TAU)
         DEALLOCATE(T)
         DEALLOCATE(Ts)

      END SUBROUTINE
