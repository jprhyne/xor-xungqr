      PROGRAM TEST_ZHT3RK
         ! Local variables
         INTEGER           I, J, K, L, N, O, II, IJ
         CHARACTER         UPLOL,UPLOR
         LOGICAL           TERMINATE
         DOUBLE PRECISION  NORM_WORK, NORM_ERR, TMP
         DOUBLE PRECISION  ALPHA, BETA
         ! Local arrays
         COMPLEX*16, ALLOCATABLE :: T(:,:), C(:,:), Ts(:,:),
     $            Cs(:,:), WORK(:,:)
         CHARACTER :: UPLOTS(2), UPLOCS(2), TRANS(2), UNITS(2)
         ! External Subroutines
         EXTERNAL ZHTRK, ZLACPY, ZLASET
         ! External Functions
         DOUBLE PRECISION  ZLANGE
         EXTERNAL          ZLANGE
         ! Intrinsic functions
         INTRINSIC SQRT,CONJG,DBLE
         ! Parameters
         COMPLEX*16 ONE, ZERO
         PARAMETER(ONE=(1.0D+0,0.0D+0), ZERO=(0.0D+0,0.0D+0))
         ! Beginning of executable statements
         TERMINATE = .FALSE.
         ! Set values for our flag vectors
         !
         ! Index with I
         UPLOTS(1) = 'U'
         UPLOTS(2) = 'L'
         ! Index with J
         UPLOCS(1) = 'U'
         UPLOCS(2) = 'L'
         ! Index with N
         TRANS(1) = 'N'
         TRANS(2) = 'C'
         ! Index with O
         UNITS(1) = 'U'
         UNITS(2) = 'N'

         ! Read in the value of K
         WRITE(*,*) "Provide a value for K to ZHT3RK"
         READ(*,*) K

         ! allocate our memory
         ALLOCATE(T(K,K))
         ALLOCATE(C(K,K))
         ALLOCATE(Ts(K,K))
         ALLOCATE(Cs(K,K))
         ALLOCATE(WORK(K,K))

         DO I = 1, 2 !UPLOTS(I): 'u','l'
         DO J = 1, 2 !UPLOCS(J): 'u','l'
         DO N = 1, 2 !TRANS(N):  'n','c'
         DO O = 1, 2 !UNITS(O):  'u','n'
            WRITE(*,*) "Flags=",UPLOTS(I), UPLOCS(J), TRANS(N),
     $                  UNITS(O)
            ! Fill our arrays with random data
            CALL RANDOM_COMPLEX_MAT(K,K,T,K)
            CALL RANDOM_COMPLEX_MAT(K,K,C,K)
            CALL RANDOM_COMPLEX_MAT(K,K,WORK,K)
            CALL RANDOM_NUMBER(ALPHA)
            CALL RANDOM_NUMBER(BETA)
            ! Copy whatever is in T and C into Ts and Cs respectively
            CALL ZLACPY('All', K, K, T, K, Ts, K)
            CALL ZLACPY('All', K, K, C, K, Cs, K)
            ! Compute \alpha*op(T)*op(T)**H + \beta*C
            CALL ZHT3RK(UPLOTS(I), UPLOCS(J), TRANS(N), UNITS(O), K,
     $            ALPHA, T, K, BETA, C, K)
            ! Check that T was not touched
            DO II = 1, K
            DO IJ = 1, K
               IF(T(II,IJ).NE.Ts(II,IJ)) THEN
                  WRITE(*,*) "T was modified at index",
     $                       "(", II,",",IJ,")"
                  TERMINATE=.TRUE.
               END IF
            END DO
            END DO
            ! If T was touched, bail. 
            IF(TERMINATE) THEN
               GOTO 10
            END IF
            ! Check that only the correct part of C was modified
            IF(UPLOCS(J).EQ.'U') THEN
               ! Check the the lower triangular part of C was not touched
               DO II = 2, K
               DO IJ = 1, II-1
                  IF(C(II,IJ).NE.Cs(II,IJ)) THEN
                     WRITE(*,*) "C was modified at index",
     $                       "(", II,",",IJ,")"
                     TERMINATE=.TRUE.
                  END IF
               END DO
               END DO
            ELSE
               ! Check the the upper triangular part of C was not touched
               DO II = 1, K-1
               DO IJ = II+1, K
                  IF(C(II,IJ).NE.Cs(II,IJ)) THEN
                     WRITE(*,*) "C was modified at index",
     $                       "(", II,",",IJ,")"
                     TERMINATE=.TRUE.
                  END IF
               END DO
               END DO
            END IF
            ! If C was touched, bail. 
            IF(TERMINATE) THEN
               GOTO 10
            END IF
            ! Determine if we performed any operation
            CALL ZLACPY('All', K, K, Cs, K, WORK, K)
            TMP = DBLE(ZERO)
            DO II = 1, K
            DO IJ = 1, K
               WORK(II,IJ) = C(II,IJ) - WORK(II,IJ)
               TMP = TMP + WORK(II,IJ)*CONJG(WORK(II,IJ))
            END DO
            END DO
            !
            IF( TMP.NE.DBLE(ZERO) ) THEN
               ! Set WORK to be 0
               CALL ZLASET('All', K, K, ZERO, ZERO, WORK, K)
               ! Copy the proper part of C into WORK
               CALL ZLACPY(UPLOCS(J), K, K, C, K, WORK, K)
               ! Set T and C to 0
               CALL ZLASET('All', K, K, ZERO, ZERO, T, K)
               CALL ZLASET('All', K, K, ZERO, ZERO, C, K)
               ! Copy the correct part of Ts and Cs into T and C respectively
               CALL ZLACPY(UPLOTS(I), K, K, Ts, K, T, K)
               CALL ZLACPY(UPLOCS(J), K, K, Cs, K, C, K)
               ! If T was supposed to be unit, then set the diagonal to ONE
               IF( UNITS(O).EQ.'U' ) THEN
                  DO II = 1, K
                     T(II,II) = ONE
                  END DO
               END IF
               ! Now, we Compute C = op(T)*op(T)**H + \beta*C using HERK
               !  since we have explicit 0s in T where needed
               CALL ZHERK(UPLOCS(J), TRANS(N), K, K, ALPHA, T, K,
     $               BETA, C, K)
               ! Now, lets compute WORK = C - WORK
               DO II = 1, K
               DO IJ = 1, K
                  WORK(II,IJ) = WORK(II,IJ) - C(II,IJ)
               END DO
               END DO
               NORM_ERR  = ZLANGE('F', K, K, WORK, K, C)
               NORM_WORK = ZLANGE('F', K, K, C, K, WORK)
               !
               NORM_ERR = NORM_ERR / NORM_WORK
               WRITE(*,*) "||actual - expected||_F / ||expected||_F = ",
     $            NORM_ERR
               IF (NORM_ERR.GT.(1.0D-10)) GOTO 10 ! Used during development. Should be uncommented in practice
            ELSE
               WRITE(*,*) "No operation was done"
            END IF
         END DO
         END DO
         END DO
         END DO

         ! Free the memory
10       DEALLOCATE(T)
         DEALLOCATE(C)
         DEALLOCATE(Ts)
         DEALLOCATE(Cs)
         DEALLOCATE(WORK)
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
