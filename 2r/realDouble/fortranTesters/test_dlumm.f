      SUBROUTINE TEST_DLUMM(N, ALPHA)
         ! Inputs
         ! Integer variables
         INTEGER           N
         ! Double variables
         DOUBLE PRECISION  ALPHA


         ! Local variables
         DOUBLE PRECISION  NORMF, TMP, TMP2

         INTEGER           I, J, K, II, IJ

         ! Local arrays
         DOUBLE PRECISION, ALLOCATABLE :: L(:,:), U(:,:), A(:,:)

         CHARACTER :: DIAGLS(2), DIAGUS(2), SIDES(2)

         ! External subroutines
         EXTERNAL DLUMM, DTRMM

         ! Parameters
         DOUBLE PRECISION  ONE, ZERO
         PARAMETER(ONE=1.0D+0, ZERO=0.0D+0)
         ! Beginning of executable statements
         DIAGLS(1) = 'U'
         DIAGLS(2) = 'N'

         DIAGUS(1) = 'U'
         DIAGUS(2) = 'N'

         SIDES(1) = 'L'
         SIDES(2) = 'R'

         ! Allocate our memory
         ALLOCATE(A(N,N))
         ALLOCATE(L(N,N))
         ALLOCATE(U(N,N))

         DO I = 1, 2 ! for each value in DIAGLS
            DO J = 1, 2 ! for each value in DIAGUS
               IF ((I.EQ.J).AND.(I.EQ.2)) THEN
                  GOTO 10 ! We do not allow both L and U to be non-unit diagonal
               END IF
               DO K = 1, 2 ! for each value in SIDES
                  ! Populate A, L, and U
                  CALL RANDOM_NUMBER(A)
                  CALL RANDOM_NUMBER(L)
                  CALL RANDOM_NUMBER(U)
                  ! Zero out the above diagonal part of L and the below diagonal part of U
                  DO II=1, N-1
                     DO IJ = II+1, N
                        L(II,IJ) = ZERO
                        U(IJ,II) = ZERO
                        ! And also copy the strictly upper triangular part of U and the strictly
                        ! lower triangular part of L into A respectively
                        A(II,IJ) = U(II,IJ)
                        A(IJ,II) = L(IJ,II)
                     END DO
                  END DO
                  ! Determine which (if either) we get the diagonal of A from
                  IF ((I.NE.J).AND.(I.EQ.1)) THEN
                     ! L is assumed unit while U has explicit values on the diagonal
                     DO II = 1, N
                        A(II,II) = U(II,II)
                     END DO
                  ELSE IF ((I.NE.J).AND.(J.EQ.1)) THEN
                     ! U is assumed unit while L has explicit values on the diagonal
                     DO II = 1, N
                        A(II,II) = L(II,II)
                     END DO
                  END IF
                  ! If neither of the above are met, then we have that I == J == 1, so we should not
                  ! be referencing the diagonal of A at all (except to set it's value).
                  ! Call our routine
                  CALL DLUMM(SIDES(K), DIAGLS(I), DIAGUS(J), N,
     $                     ALPHA, A, N)

                  ! Now, we want to compute the expected product, which is just a TRMM that will
                  ! be stored in U. So we will set the diagonal of U if needed
                  IF (J.EQ.1) THEN
                     DO II = 1, N
                        U(II,II) = ONE
                     END DO
                  END IF
                  CALL DTRMM(SIDES(K), 'Lower', 'No Transpose',
     $                  DIAGLS(I), N, N, ALPHA, L, N, U, N)

                  ! Now, we will compute the relative forward error in the frobenius norm
                  TMP2 = ZERO
                  NORMF= ZERO
                  DO II = 1, N
                     DO IJ = 1, N
                        TMP = A(II,IJ) - U(II,IJ)
                        NORMF = NORMF + TMP*TMP
                        TMP2 = TMP2 + U(II,IJ)*U(II,IJ)
                     END DO
                  END DO
                  IF (TMP2.NE.ZERO) THEN
                     NORMF = NORMF / TMP2
                  END IF
                  NORMF = SQRT(NORMF)
                  WRITE(*,*) "Parameters to DLUMM"
                  WRITE(*,*) "SIDE=",SIDES(K)," DIAGL=",DIAGLS(I),
     $               " DIAGU=",DIAGUS(J)
                  ! Now, print the error
                  WRITE(*,*) NORMF
               END DO
10             CONTINUE
            END DO
         END DO
         ! Free our memory
         DEALLOCATE(A)
         DEALLOCATE(L)
         DEALLOCATE(U)
      END SUBROUTINE
