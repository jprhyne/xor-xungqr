      SUBROUTINE TESTDTRMMOOP(M, N, ALPHA, BETA)
         ! Arguments
         INTEGER  M, N
         DOUBLE PRECISION  ALPHA, BETA

         ! Local variables
         INTEGER     I,J,K,L,MAXMN, II, IJ, O
         LOGICAL   TERMINATE
         DOUBLE PRECISION  NORM_F, TMP
         ! Local arrays
         DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:), C(:,:),
     $            As(:,:), Bs(:,:), Cs(:,:), WORK(:,:)
         CHARACTER :: SIDES(2), TRANSAS(3), TRANSBS(3), DIAGS(2),
     $            UPLOS(2)
         ! External Subroutines
         EXTERNAL DTRMMOOP, DLACPY
         ! External Functions
         ! Intrinsic functions
         INTRINSIC MAX
         ! Parameters
         DOUBLE PRECISION ZERO
         PARAMETER(ZERO=0.0D+0)
         ! Beginning of executable statements
         TERMINATE = .FALSE.

         SIDES(1) = 'L'
         SIDES(2) = 'R'

         TRANSAS(1) = 'N'
         TRANSAS(2) = 'T'
         TRANSAS(3) = 'C'

         TRANSBS(1) = 'N'
         TRANSBS(2) = 'T'
         TRANSBS(3) = 'C'

         DIAGS(1) = 'N'
         DIAGS(2) = 'U'

         UPLOS(1) = 'U'
         UPLOS(2) = 'L'

         MAXMN = M
         IF (N.GT.M) MAXMN = N

         ! Allocate our memory
         ALLOCATE(A(MAXMN,MAXMN))
         ALLOCATE(As(MAXMN,MAXMN))
         ALLOCATE(B(MAXMN,MAXMN))
         ALLOCATE(Bs(MAXMN,MAXMN))
         ALLOCATE(C(MAXMN,MAXMN))
         ALLOCATE(Cs(MAXMN,MAXMN))
         ALLOCATE(WORK(MAXMN,MAXMN))

         DO I = 1, 2 ! For each side
            DO J = 1, 3 ! For each transposing of A
               DO K = 1, 3 ! For each transposing of B
                  DO L = 1, 2 ! For each of B being unit or non-unit
                     DO O = 1, 2 ! For each of B being upper or lower
                        ! Fill our matrices up with random numbers
                        CALL RANDOM_NUMBER(A)
                        CALL RANDOM_NUMBER(B)
                        CALL RANDOM_NUMBER(C)
                        CALL RANDOM_NUMBER(WORK)

                        ! Store A,B,C in As,Bs,Cs respectively
                        CALL DLACPY('All', MAXMN, MAXMN, A, MAXMN, As,
     $                     MAXMN)
                        CALL DLACPY('All', MAXMN, MAXMN, B, MAXMN, Bs,
     $                     MAXMN)
                        CALL DLACPY('All', MAXMN, MAXMN, C, MAXMN, Cs,
     $                     MAXMN)

                        ! Compute our desired operation
                        CALL DTRMMOOP(SIDES(I), UPLOS(O), TRANSAS(J),
     $                     TRANSBS(K), DIAGS(L), M, N, ALPHA, A, MAXMN,
     $                     B, MAXMN, BETA, C, MAXMN)

                        ! Make sure that A and B were not modified
                        DO II = 1, MAXMN
                           DO IJ = 1, MAXMN
                              IF (A(II,IJ).NE.As(II,IJ)) THEN
                                 WRITE(*,*) "A was modified at index",
     $                              "(",II,",",IJ,")"
                                 TERMINATE = .TRUE.
                              END IF
                              IF (B(II,IJ).NE.Bs(II,IJ)) THEN
                                 WRITE(*,*) "B was modified at index",
     $                              "(",II,",",IJ,")"
                                 TERMINATE = .TRUE.
                              END IF
                           END DO
                        END DO
                        IF (TERMINATE) THEN
                           GOTO 10 ! Make sure we free our memory
                        END IF
                        ! Determine if we did the correct operations
                        IF (TRANSAS(J).EQ.'N') THEN
                           CALL DLACPY('ALL', MAXMN, MAXMN, As, MAXMN,
     $                        WORK, MAXMN)
                        ELSE
                           DO II = 1, MAXMN
                              DO IJ = 1, MAXMN
                                 WORK(II,IJ) = As(IJ,II)
                              END DO
                           END DO
                        END IF
                        CALL DTRMM(SIDES(I), UPLOS(O), TRANSBS(K),
     $                     DIAGS(L), M, N, ALPHA, B, MAXMN, WORK, MAXMN)
                        ! Now add \beta C to WORK
                        DO II = 1, M
                           DO IJ = 1, N
                              WORK(II,IJ) = WORK(II,IJ) + BETA*Cs(II,IJ)
                           END DO
                        END DO
                        ! Print out the error to console for visual inspection
                        NORM_F = 0.0
                        DO II = 1, M
                           DO IJ = 1, N
                              TMP = WORK(II,IJ) - C(II,IJ)
                              NORM_F = NORM_F + TMP * TMP
                           END DO
                        END DO
                        NORM_F = SQRT(NORM_F)
                        TMP = 0.0
                        DO II = 1, M
                           DO IJ = 1, N
                              TMP = TMP + WORK(II,IJ) * WORK(II,IJ)
                           END DO
                        END DO
                        TMP = SQRT(TMP)
                        IF (TMP.NE.ZERO) THEN
                           NORM_F = NORM_F / TMP
                        END IF
                        ! Print out the flags used to allow for repeatability
                        WRITE(*,*) "Parameters to DTRMMOOP"
                        WRITE(*,*) "Side=",SIDES(I)," UPLO=",UPLOS(O),
     $                     " TRANSA=",TRANSAS(J)," TRANSB=", TRANSBS(K),
     $                     " DIAG=", DIAGS(L)
                        ! Print the error out
                        WRITE(*,*) "Forward error: ", NORM_F
                     END DO
                  END DO
               END DO
            END DO
         END DO

         ! Free our memory
10       DEALLOCATE(A)
         DEALLOCATE(As)
         DEALLOCATE(B)
         DEALLOCATE(Bs)
         DEALLOCATE(C)
         DEALLOCATE(Cs)
         DEALLOCATE(WORK)

      END SUBROUTINE
