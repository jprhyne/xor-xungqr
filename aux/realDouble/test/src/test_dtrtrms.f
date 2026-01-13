      PROGRAM TEST_DTRTRMS
      INTEGER     I,J,K,L,M, N, INFO, II, IJ
      DOUBLE PRECISION  NORM_ERR, NORM_X, ONE, ZERO
      CHARACTER         SIDET, UPLO, TRANS, DIAGT, DIAGV
      LOGICAL           TUNIT, VUNIT, TERMINATE, TRANSV, TUPPER, VUPPER

      DOUBLE PRECISION, ALLOCATABLE :: T(:,:), T_INV(:,:), WORKMAT(:,:),
     $      V(:,:), Vs(:,:), Ts(:,:)
      CHARACTER :: SIDETS(2), UPLOS(2), TRANSS(3), DIAGTS(2), DIAGVS(2)

      DOUBLE PRECISION  DLANGE
      EXTERNAL          DLANGE

      EXTERNAL          DTRTRMS, DTRTRI, DTRMM, DLACPY

      ! Read in our dimension
      READ(*,*) N

      ONE = 1.0D+0
      ZERO = 0.0D+0

      ! Allocate our memory
      ALLOCATE(T(N,N))
      ALLOCATE(V(N,N))
      ALLOCATE(Ts(N,N))
      ALLOCATE(Vs(N,N))
      ALLOCATE(T_INV(N,N))
      ALLOCATE(WORKMAT(N,N))

      SIDETS(1) = "L"
      SIDETS(2) = "R"

      UPLOS(1) = "U"
      UPLOS(2) = "L"

      TRANSS(1) = "N"
      TRANSS(2) = "T"
      TRANSS(3) = "C"

      DIAGTS(1) = "U"
      DIAGTS(2) = "N"

      DIAGVS(1) = "U"
      DIAGVS(2) = "N"

      DO I = 1,2 ! SIDETS = ( "L", "R" )
      DO J = 1,2 ! UPLOS  = ( "U", "L" )
      DO K = 1,3 ! TRANSS = ( "N", "T", "C" )
      DO L = 1,2 ! DIAGTS = ( "U", "N" )
      DO M = 1,2 ! DIAGVS = ( "U", "N" )
         SIDET = SIDETS(I)
         UPLO = UPLOS(J)
         TRANS = TRANSS(K)
         DIAGT = DIAGTS(L)
         DIAGV = DIAGVS(M)
         WRITE(*,*) "Running Flags: ", SIDET, UPLO, TRANS, DIAGT, DIAGV

         TUNIT = DIAGT.EQ.'U'
         VUNIT = DIAGV.EQ.'U'

         TUPPER = UPLO.EQ.'U'
         TRANSV = TRANS.EQ.'T'.OR.TRANS.EQ.'C'

         ! Determing if V is stored as upper or lower
         VUPPER = (TUPPER.AND.(.NOT.TRANSV)).OR.
     $               ((.NOT.TUPPER).AND.TRANSV)
         TERMINATE = .FALSE.
         ! Generate our T and V matrices
         CALL RANDOM_NUMBER(T)
         CALL RANDOM_NUMBER(V)
         ! Add one to the diagonal to move away from singularity
         DO II = 1, N
            V(II,II) = V(II,II) + ONE
            T(II,II) = T(II,II) + ONE
         END DO
         ! Copy V into Vs
         CALL DLACPY('All', N, N, V, N, Vs, N)
         ! Copy the correct part
         CALL DLACPY(UPLO, N, N, T, N, T_INV, N)
         ! Compute X such that TX = op(V) or XT = op(V)
         CALL DTRTRMS(SIDET, UPLO, TRANS, DIAGT, DIAGV, N, T, N, V, N)
         ! Ensure that V was not modified. If it was, bail immediately
         DO II = 1, N
            DO IJ = 1, N
               IF (V(II,IJ) .NE. Vs(II,IJ)) THEN
                  TERMINATE = .TRUE.
                  WRITE(*,*) "V(",II,",",IJ,") was modified"
               END IF
            END DO
         END DO
         IF (TERMINATE) GOTO 10
         ! Compute T^-1
         CALL DTRTRI(UPLO, DIAGT, N, T_INV, N, INFO)
         ! Compute WORKMAT = T^{-1}op(V) or op(V)T^{-1}
         ! Set V to 0
         CALL DLASET('All', N, N, ZERO, ZERO, V, N)
         ! Copy the correct part of Vs into V
         IF (VUPPER) THEN
            CALL DLACPY('U', N, N, Vs, N, V, N)
         ELSE
            CALL DLACPY('L', N, N, Vs, N, V, N)
         END IF

         IF (VUNIT) THEN
            DO II = 1, N
               V(II,II) = ONE
            END DO
         END IF
         ! copy op(V) into workmat
         IF (TRANSV) THEN
            DO II = 1, N
            DO IJ = 1, N
               WORKMAT(IJ,II) = V(II,IJ)
            END DO
            END DO
         ELSE
            CALL DLACPY('All', N, N, V, N, WORKMAT, N)
         END IF

         ! Compute WORKMAT = WORKMAT*T^{-1} or T^{-1}*WORKMAT
         CALL DTRMM(SIDET, UPLO, 'N', DIAGT, N, N, ONE, T_INV, N,
     $         WORKMAT, N)
         ! Compute NORM_X = ||WORKMAT||_F
         NORM_X = DLANGE('F', N, N, WORKMAT, N, WORKMAT)
         ! Compute WORKMAT = WORKMAT - X = WORKMAT - T
         IF (TUPPER) THEN
            DO IJ = 1, N
            DO II = 1, IJ
               WORKMAT(II,IJ) = WORKMAT(II,IJ) - T(II,IJ)
            END DO
            END DO
         ELSE
            DO IJ = 1, N
            DO II = IJ, N
               WORKMAT(II,IJ) = WORKMAT(II,IJ) - T(II,IJ)
            END DO
            END DO
         END IF
         ! Compute NORM_ERR = ||WORKMAT||_F/ NORM_X
         NORM_ERR = DLANGE('F', N, N, WORKMAT, N, WORKMAT) / NORM_X
         ! Print NORM_ERR = ||X_MULT - X_MINE||_F / ||X_MINE||_F
         WRITE(*,*) NORM_ERR
      END DO
      END DO
      END DO
      END DO
      END DO

10    DEALLOCATE(T)
      DEALLOCATE(V)
      DEALLOCATE(Ts)
      DEALLOCATE(Vs)
      DEALLOCATE(T_INV)
      DEALLOCATE(WORKMAT)
      END PROGRAM
