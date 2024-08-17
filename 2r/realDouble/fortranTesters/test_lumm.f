      SUBROUTINE TEST_LUMM(N, LDL, LDU)
         INTEGER N, LDL, LDU

         ! Local variables
         DOUBLE PRECISION  NORMF, TMP, ALPHA

         INTEGER           I, J, LDQ

         ! Local arrays
         DOUBLE PRECISION, ALLOCATABLE :: L(:,:), U(:,:)
         DOUBLE PRECISION, ALLOCATABLE :: Ls(:,:), Us(:,:)
         DOUBLE PRECISION, ALLOCATABLE :: A(:,:), Q(:,:)

         ! External subroutines
         EXTERNAL LUMM, DLACPY, DTRMM

         ! Parameters
         DOUBLE PRECISION  ONE, ZERO
         PARAMETER(ONE=1.0D+0, ZERO=0.0D+0)
         LDQ = 2*MAX(LDL,LDU)
         ! Allocate memory
         ALLOCATE(L(LDL,N))
         ALLOCATE(U(LDU,N))
         ALLOCATE(Ls(N,N))
         ALLOCATE(Us(N,N))
         ALLOCATE(A(N,N))
         ALLOCATE(Q(LDQ,N))

         ! Generate L and U as random matrices
         CALL RANDOM_NUMBER(L)
         CALL RANDOM_NUMBER(U)

         DO I = 1, N
            DO J = 1, N
               IF(I.LT.J) THEN
                  L(I,J) = 0
               ELSE IF (I.GT.J) THEN
                  U(I,J) = 0
               END IF
            END DO
         END DO

         CALL DLACPY('Lower', N, N, L, LDL, Q, LDQ)
         CALL DLACPY('Upper', N, N, U, LDU, Q, LDQ)

         CALL DLACPY('All', N, N, L, LDL, Ls, N)
         CALL DLACPY('All', N, N, U, LDU, Us, N)
         ALPHA = 1.0D+0


         CALL LUMM(N, ALPHA, Q, LDQ)

         ! Using Ls and Us compute Ls * Us using TRMM
         ! A = Ls, B = Us
         ! First, since we are treating Us as the 'non-triangular' matrix in 
         ! dtrmm, we need to first enforce it being upper triangular.
         DO I = 2, N
            DO J = 1, I - 1
               Us(I,J) = 0
            END DO
         END DO
         CALL DTRMM('Left', 'Lower', 'No', 'Unit', N, N, ALPHA,
     $      Ls, N, Us, N)

         ! Compute Us - Q and find the norm thereof
         NORMF=0
         DO I = 1, N
            DO J = 1, N
               TMP = Us(I,J) - Q(I,J)
               NORMF = NORMF + TMP * TMP
            END DO
         END DO

         WRITE(*,*) "ALPHA = ", ALPHA, ": ||Us - Q||_F = ", NORMF

         ! Try with alpha = -1, as this is the functionality we need
         ALPHA = -1.0D+0
         ! Generate L and U as random matrices
         CALL RANDOM_NUMBER(L)
         CALL RANDOM_NUMBER(U)

         DO I = 1, N
            DO J = 1, N
               IF(I.LT.J) THEN
                  L(I,J) = 0
               ELSE IF (I.GT.J) THEN
                  U(I,J) = 0
               END IF
            END DO
         END DO

         CALL DLACPY('Lower', N, N, L, LDL, Q, LDQ)
         CALL DLACPY('Upper', N, N, U, LDU, Q, LDQ)

         CALL DLACPY('All', N, N, L, LDL, Ls, N)
         CALL DLACPY('All', N, N, U, LDU, Us, N)

         CALL LUMM(N, ALPHA, Q, LDQ)

         ! Using Ls and Us compute Ls * Us using TRMM
         ! A = Ls, B = Us
         ! First, since we are treating Us as the 'non-triangular' matrix in 
         ! dtrmm, we need to first enforce it being upper triangular.
         DO I = 2, N
            DO J = 1, I - 1
               Us(I,J) = 0
            END DO
         END DO
         CALL DTRMM('Left', 'Lower', 'No', 'Unit', N, N, ALPHA,
     $      Ls, N, Us, N)

         ! Compute Us - Q and find the norm thereof
         NORMF=0
         DO I = 1, N
            DO J = 1, N
               TMP = Us(I,J) - Q(I,J)
               NORMF = NORMF + TMP * TMP
            END DO
         END DO

         WRITE(*,*) "ALPHA = ", ALPHA, ": ||Us - Q||_F = ", NORMF

         DEALLOCATE(L)
         DEALLOCATE(Ls)
         DEALLOCATE(U)
         DEALLOCATE(Us)
         DEALLOCATE(A)
         DEALLOCATE(Q)
      END SUBROUTINE
