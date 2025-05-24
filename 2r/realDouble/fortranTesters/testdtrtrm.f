      SUBROUTINE TESTDTRTRM(N, ALPHA)
         ! Arguments
         INTEGER  N
         DOUBLE PRECISION  ALPHA

         ! Local variables
         INTEGER           I, J, K, L, O, II, IJ
         CHARACTER         UPLOV
         LOGICAL           TERMINATE, VLOWER
         DOUBLE PRECISION  NORM_WORK, NORM_ERR, TMP
         ! Local arrays
         DOUBLE PRECISION, ALLOCATABLE :: T(:,:), V(:,:), Ts(:,:),
     $            Vs(:,:), WORK(:,:)
         CHARACTER :: SIDES(2), UPLOS(2), TRANSVS(3), DIAGTS(2),
     $            DIAGVS(2)
         ! External Subroutines
         EXTERNAL DTRTRM, DLACPY
         ! External Functions
         DOUBLE PRECISION  DLANGE
         EXTERNAL          DLANGE
         ! Intrinsic functions
         INTRINSIC SQRT
         ! Parameters
         DOUBLE PRECISION ONE, ZERO
         PARAMETER(ONE=1.0D+0, ZERO=0.0D+0)
         ! Beginning of executable statements
         TERMINATE = .FALSE.
         ! Set our arrays
         SIDES(1) = 'L'
         SIDES(2) = 'R'

         UPLOS(1) = 'U'
         UPLOS(2) = 'L'

         TRANSVS(1) = 'T'
         TRANSVS(2) = 'C'
         TRANSVS(3) = 'N'

         DIAGTS(1) = 'U'
         DIAGTS(2) = 'N'

         DIAGVS(1) = 'U'
         DIAGVS(2) = 'N'
         ! Allocate our memory
         ALLOCATE(T(N,N))
         ALLOCATE(Ts(N,N))
         ALLOCATE(V(N,N))
         ALLOCATE(Vs(N,N))
         ALLOCATE(WORK(N,N))

         DO I = 1, 2 ! For each side op(V) is on
            DO J = 1, 2 ! For each of T being upper and lower triangular
               DO K = 1, 3 ! For each of V being Transposed
                  DO L = 1, 2 ! For each of T being unit or not
                     DO O = 1, 2 ! For each of V being unit or not
                        ! Fill our matrices up with random numbers
                        CALL RANDOM_NUMBER(T)
                        CALL RANDOM_NUMBER(V)
                        ! Copy V and T into Vs and Ts respectively
                        CALL DLACPY('All', N, N, V, N, Vs, N)
                        CALL DLACPY('All', N, N, T, N, Ts, N)
                        ! Compute our matrix product
                        CALL DTRTRM(SIDES(I), UPLOS(J), TRANSVS(K),
     $                        DIAGTS(L), DIAGVS(O), N, ALPHA, T, N,
     $                        V, N)
                        ! Check to make sure that V was not changed
                        DO II = 1, N
                           DO IJ = 1, N
                              IF(V(II,IJ).NE.Vs(II,IJ)) THEN
                                 WRITE(*,*) "V was modified at index",
     $                              "(", II,",",IJ,")"
                                 TERMINATE=.TRUE.
                              END IF
                           END DO
                        END DO
                        ! Check to make sure only the part of T we are
                        ! considering is touched
                        !IF (UPLOS(J).EQ.'U') THEN
                           ! Check that we did not touch the lower part
                        IF (TERMINATE) THEN
                           GOTO 10 ! Free our memory then exit!
                        END IF
                        ! Set work to be the 0 matrix explicitly
                        CALL DLASET('All', N, N, ZERO, ZERO, WORK, N)
                        ! Now, make sure that we did the correct operations
                        ! Copy the correct triangular part of Ts WORK
                        CALL DLACPY(UPLOS(J), N, N, Ts, N, WORK, N)
                        ! If T was assumed unit triangular, set the
                        ! diagonal to 1 explicitly
                        IF (DIAGTS(L).EQ.'U') THEN
                           DO II = 1, N
                              WORK(II,II) = ONE
                           END DO
                        END IF
                        ! The shape of V will be dependent on UPLO and TRANSV
                        VLOWER = UPLOS(J).EQ.'L'.AND.TRANSVS(K).EQ.'N'
                        VLOWER = VLOWER.OR.(UPLOS(J).EQ.'U'.AND.
     $                              TRANSVS(K).EQ.'C')
                        VLOWER = VLOWER.OR.(UPLOS(J).EQ.'U'.AND.
     $                              TRANSVS(K).EQ.'T')
                        IF(VLOWER) THEN
                           UPLOV = 'L'
                        ELSE
                           UPLOV = 'U'
                        END IF
                        ! Now, compute the correct operation
                        CALL DTRMM(SIDES(I), UPLOV, TRANSVS(K),
     $                           DIAGVS(O), N, N, ALPHA, Vs, N, WORK, N)
                        ! Compute the relative forward error
                        NORM_WORK = DLANGE('F', N, N, WORK, N, WORK, N)
                        ! Set the part of T that was not supposed to be
                        ! touched to 0 to make our error computation easier
                        IF (UPLOS(J).EQ.'L') THEN
                           DO II = 1, N-1
                              DO IJ = II+1, N
                                 T(II,IJ) = ZERO
                              END DO
                           END DO
                        ELSE
                           DO II = 2, N
                              DO IJ = 1, II-1
                                 T(II,IJ) = ZERO
                              END DO
                           END DO
                        END IF
                        NORM_ERR = ZERO
                        DO II = 1, N
                           DO IJ = 1, N
                              TMP = T(II, IJ) - WORK(II, IJ)
                              NORM_ERR = NORM_ERR + TMP*TMP
                           END DO
                        END DO
                        ! Don't divide by 0
                        IF (NORM_WORK.NE.ZERO) THEN
                           NORM_ERR = NORM_ERR/NORM_WORK
                        END IF
                        NORM_ERR = SQRT(NORM_ERR)
                        ! Print out the flags used to allow for repeatability
                        WRITE(*,*) "Parameters to DTRTRM"
                        WRITE(*,*) "Side=",SIDES(I)," UPLO=",UPLOS(J),
     $                     " TRANSV=",TRANSVS(K)," DIAGT=", DIAGTS(L),
     $                     " DIAGV=", DIAGVS(O)
                        ! Print the error out
                        WRITE(*,*) "Forward error: ", NORM_ERR
                     END DO
                  END DO
               END DO
            END DO
         END DO

         ! Free our memory
10       DEALLOCATE(T)
         DEALLOCATE(Ts)
         DEALLOCATE(V)
         DEALLOCATE(Vs)
         DEALLOCATE(WORK)

      END SUBROUTINE
