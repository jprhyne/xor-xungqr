      PROGRAM TEST_ZTRMMOOP_LVL2
         ! Scalars
         INTEGER           I,J,K,L,M,N,O,MN,II,IJ
         CHARACTER         VERB_CHAR
         LOGICAL           TERMINATE,VERBOSE
         DOUBLE PRECISION  NORM_ERR, NORM_C
         COMPLEX*16        ALPHA, BETA, TMP
         ! Arrays
         COMPLEX*16, ALLOCATABLE :: A(:,:), As(:,:), B(:,:),
     $      Bs(:,:), C(:,:), Cs(:,:)
         CHARACTER :: SIDES(2), UPLOS(2), TRANSAS(3), TRANSBS(3),
     $      DIAGS(2)
         ! External subroutines
         EXTERNAL          ZTRMMOOP, ZTRMMOOP_LVL2, ZLACPY
         ! External Functions
         DOUBLE PRECISION  ZLANGE
         EXTERNAL          ZLANGE
         ! Intrinsic Functions
         INTRINSIC         MAX
         ! Parameters
         DOUBLE PRECISION  DONE, DZERO
         COMPLEX*16        ONE, ZERO
         PARAMETER(DONE=1.0D+0, DZERO=0.0D+0)
         PARAMETER(ONE=(1.0D+0,0.0D+0), ZERO=(0.0D+0,0.0D+0))
         ! Beginning of executable statements
         TERMINATE = .FALSE.
         VERBOSE = .FALSE.
         ! Set our flag vectors
         ! Index with I
         SIDES(1) = 'L'
         SIDES(2) = 'R'
         ! Index with J
         UPLOS(1) = 'L'
         UPLOS(2) = 'U'
         ! Index with K
         TRANSAS(1) = 'N'
         TRANSAS(2) = 'T'
         TRANSAS(3) = 'C'
         ! Index with L
         TRANSBS(1) = 'N'
         TRANSBS(2) = 'T'
         TRANSBS(3) = 'C'
         ! Index with O
         DIAGS(1) = 'U'
         DIAGS(2) = 'N'

         ! Read in the values for M and N
         WRITE(*,*) "Provide a value for M and N to ZTRMMOOP_LVL2"
         READ(*,*) M, N
         WRITE(*,*) "You provided (M,N) = (",M,",",N,")"
         WRITE(*,*) "Do you want 'no op' information? 1=yes, otherwise",
     $      " no"
         READ(*,*) VERB_CHAR
         VERBOSE = VERB_CHAR.EQ.'1'
         IF (VERBOSE) THEN
            WRITE(*,*) "no op information will be printed"
         ELSE
            WRITE(*,*) "no op information will NOT be printed"
         END IF

         ! Allocate our memory. Determine our max dimension for B
         MN = MAX(M,N)
         ! If MIN(M,N) < 0 then we return immediately
         IF (MIN(M,N).LT.0) THEN
            WRITE(*,*) "You gave a negative input value dummy-head"
            RETURN
         END IF
         ALLOCATE(A(MN,MN))
         ALLOCATE(As(MN,MN))
         ALLOCATE(B(MN,MN))
         ALLOCATE(Bs(MN,MN))
         ALLOCATE(C(M,N))
         ALLOCATE(Cs(M,N))
         !RLNNU
         DO I = 1, 2 ! SIDES   = ['L', 'R']
         DO J = 1, 2 ! UPLOS   = ['L', 'U']
         DO K = 1, 3 ! TRANSAS = ['N', 'T', 'C']
         DO L = 1, 3 ! TRANSBS = ['N', 'T', 'C']
         DO O = 1, 2 ! DIAGS   = ['U', 'N']
            IF (K.NE.3.AND.L.NE.3) THEN
               CONTINUE
            END IF
            ! Generate our matrices, A, B, and C as random
            CALL RANDOM_COMPLEX_MAT(MN, MN, A, MN)
            CALL RANDOM_COMPLEX_MAT(MN, MN, B, MN)
            CALL RANDOM_COMPLEX_MAT(M, N, C, M)
            ! Store these matrices into As,Bs, and Cs respectively
            CALL ZLACPY('All', MN, MN, A, MN, As, MN)
            CALL ZLACPY('All', MN, MN, B, MN, Bs, MN)
            CALL ZLACPY('All', M, N, C, M, Cs, M)
            ! Generate our scalars as random numbers
            CALL RANDOM_COMPLEX_SCALAR(ALPHA)
            CALL RANDOM_COMPLEX_SCALAR(BETA)
            ! Compute either C = alpha*op(A)*op(B) + beta*C
            ! or
            ! C = alpha*op(B)*op(A) + beta*C
            ! depending on the value of SIDES(I) using our level 2
            ! implementation
            CALL ZTRMMOOP_LVL2(SIDES(I), UPLOS(J), TRANSAS(K),
     $            TRANSBS(L), DIAGS(O), M, N, ALPHA, A, MN, B, MN, BETA,
     $            C, M)
            ! Ensure that A and B were not changed
            DO II = 1, MN
            DO IJ = 1, MN
               IF(A(II,IJ).NE.As(II,IJ)) THEN
                  WRITE(*,*) "A was modified at index",
     $               "(",II,",",IJ,")"
                  TERMINATE = .TRUE.
               END IF
            END DO
            END DO
            ! If A was touched, bail.
            IF(TERMINATE) THEN
               GOTO 10
            END IF
            DO II = 1, MN
            DO IJ = 1, MN
               IF(B(II,IJ).NE.Bs(II,IJ)) THEN
                  WRITE(*,*) "B was modified at index",
     $               "(",II,",",IJ,")"
                  TERMINATE = .TRUE.
               END IF
            END DO
            END DO
            ! If B was touched, bail.
            IF(TERMINATE) THEN
               GOTO 10
            END IF
            ! To make development easier, we check if C was actually modified
            ! and only do our checks if C was modified.
            TMP = ZERO
            NORM_ERR = DZERO
            DO II = 1, M
            DO IJ = 1, N
               TMP = C(II,IJ) - Cs(II,IJ)
               NORM_ERR = NORM_ERR + DBLE(DCONJG(TMP)*TMP)
            END DO
            END DO
            ! NORM_ERR is ZERO only when
            ! ALPHA = ZERO and BETA = ONE or
            ! no operation was done. We don't treat the first as a real option
            ! as ALPHA and BETA were randomly made
            IF (NORM_ERR.NE.ZERO) THEN
               ! Print our flags that we ran
               WRITE(*,*) "Flags=", SIDES(I), UPLOS(J), TRANSAS(K),
     $            TRANSBS(L), DIAGS(O)
               CALL ZTRMMOOP(SIDES(I), UPLOS(J), TRANSAS(K), TRANSBS(L),
     $               DIAGS(O), M, N, ALPHA, A, MN, B, MN, BETA, Cs, M)
               ! Compute the norm of Cs as this is the "true" value
               NORM_C = ZLANGE('F', M, N, Cs, M, Cs)
               ! Compute C = C - Cs
               DO II = 1, M
               DO IJ = 1, N
                  C(II,IJ) = C(II,IJ) - Cs(II,IJ)
               END DO
               END DO
               NORM_ERR = ZLANGE('F', M, N, C, M, C) / NORM_C
               WRITE(*,*) "||actual - expected||_F / ||expected||_F = ",
     $            NORM_ERR
            ELSE IF (VERBOSE) THEN! Nothing was done
               ! We print that no operation was done if the user
               ! requested
               WRITE(*,*) "Flags=", SIDES(I), UPLOS(J), TRANSAS(K),
     $            TRANSBS(L), DIAGS(O)
               WRITE(*,*) "No operation was done"
            END IF
         END DO ! End I
         END DO ! End J
         END DO ! End K
         END DO ! End L
         END DO ! End O
10       DEALLOCATE(A)
         DEALLOCATE(As)
         DEALLOCATE(B)
         DEALLOCATE(Bs)
         DEALLOCATE(C)
         DEALLOCATE(Cs)
      END PROGRAM

      SUBROUTINE RANDOM_COMPLEX_MAT(M, N, A, LDA)
         INTEGER M,N,LDA
         COMPLEX*16 A(LDA,*)

         DOUBLE PRECISION RP,CP
         INTEGER  I,J

         DO I = 1, M
         DO J = 1, N
            CALL RANDOM_NUMBER(RP)
            CALL RANDOM_NUMBER(CP)
            A(I,J) = CMPLX(RP,CP)
         END DO
         END DO
      END SUBROUTINE

      SUBROUTINE RANDOM_COMPLEX_SCALAR(A)
         COMPLEX*16 A
         DOUBLE PRECISION RP,CP
         CALL RANDOM_NUMBER(RP)
         CALL RANDOM_NUMBER(CP)
         A = CMPLX(RP,CP)
      END SUBROUTINE
