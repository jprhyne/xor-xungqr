      PROGRAM TEST_ZTRMVOOP
         ! Scalars
         INTEGER           I,J,K,L,N,INCX,II,IJ
         LOGICAL           TERMINATE
         DOUBLE PRECISION  NORM_ERR
         COMPLEX*16        ALPHA, BETA
         ! Arrays
         COMPLEX*16      , ALLOCATABLE :: A(:,:),As(:,:),X(:,:),Xs(:,:),
     $         Y(:),Ys(:),WORK(:,:)
c        CHARACTER ::      UPLOS(2), TRANS(3), DIAGS(2), DIR(4)
         CHARACTER ::      UPLOS(2), TRANS(3), DIAGS(2), DIR(2)
         ! External Subroutines
         EXTERNAL          ZTRMV, ZTRMVOOP
         ! External Functions
         DOUBLE PRECISION  DLANGE,DNRM2
         EXTERNAL          DLANGE,DNRM2
         ! Parameters
         COMPLEX*16        ONE, ZERO, NEG_ONE
         PARAMETER(ONE=(1.0D+0,0.0D+0), ZERO=(0.0D+0,0.0D+0),
     $      NEG_ONE=(-1.0D+0,0.0D+0))
         ! Beginning of executable statements
         TERMINATE = .FALSE.
         ! Set the values for our flag vectors
         !
         ! Index with I
         UPLOS(1) = 'U'
         UPLOS(2) = 'L'
         !
         ! Index with J
         TRANS(1) = 'N'
         TRANS(2) = 'T'
         TRANS(3) = 'C'
         !
         ! Index with K
         DIAGS(1) = 'U'
         DIAGS(2) = 'N'
         !
         ! Index with L
         DIR(1) = 'C' ! Column
         DIR(2) = 'R' ! Row
c        DIR(3) = 'NC' ! Negative column
c        DIR(4) = 'NR' ! Negative row

         ! Read in a value of N
         WRITE(*,*) "Provide a value for N to ZTRMVOOP"
         READ(*,*) N

         ! Allocate our memory
         ALLOCATE(A(N,N))
         ALLOCATE(As(N,N))
         ALLOCATE(X(N,N))
         ALLOCATE(Xs(N,N))
         ALLOCATE(Y(N))
         ALLOCATE(Ys(N))
         ALLOCATE(WORK(N,N))
         ! Iterate over all our possible flags
         DO I = 1, 2 ! UPLOS(I): 'U','L'
         DO J = 1, 3 ! TRANS(J): 'N','T','C'
         DO K = 1, 2 ! DIAGS(K): 'U','N'
         DO L = 1, 2 ! DIR(L): 'C','R'
            IF (DIR(L).EQ.'C') THEN
               INCX = 1
            ELSE IF (DIR(L).EQ.'R') THEN
               INCX = N
            ELSE IF (DIR(L).EQ.'NC') THEN
               !INCX = -1
            ELSE
               !INCX = -N
            END IF
            WRITE(*,*) "Flags=",UPLOS(I),TRANS(J),DIAGS(K),DIR(L)
            ! Fill our arrays with random data
            CALL RANDOM_COMPLEX_MAT(N, N, A, N)
            CALL RANDOM_COMPLEX_MAT(N, N, X, N)
            CALL RANDOM_COMPLEX_MAT(N, 1, Y, 1)
            CALL RANDOM_COMPLEX_SCALAR(ALPHA)
            CALL RANDOM_COMPLEX_SCALAR(BETA)
            ! Copy whatever is in A and Y into As and Ys respectively
            CALL ZLACPY('All', N, N, A, N, As, N)
            CALL ZLACPY('All', N, N, X, N, Xs, N)
            CALL ZCOPY(N, Y, 1, Ys, 1)
            ! Compute \alpha*op(A)x + \beta*y
            CALL ZTRMVOOP(UPLOS(I), TRANS(J), DIAGS(K), N, ALPHA, A, N,
     $            X, INCX, BETA, Y, 1)
            ! Ensure that A and X are not modified
            DO II = 1, N
            DO IJ = 1, N
               IF (A(II,IJ).NE.As(II,IJ)) THEN
                  WRITE(*,*) "A was modified at index (I,J) = (",I,",",
     $               J,")"
                  TERMINATE = .TRUE.
               END IF
            END DO
            END DO
            IF(TERMINATE) GOTO 10
            DO II = 1, N
            DO IJ = 1, N
               IF (X(II,IJ).NE.Xs(II,IJ)) THEN
                  WRITE(*,*) "X was modified at index (I,J) = (",I,",",
     $               J,")"
                  TERMINATE = .TRUE.
               END IF
            END DO
            END DO
            IF(TERMINATE) GOTO 10
            ! Now, we test the functionality by using dtrmv
            ! work = x
            CALL ZCOPY(N, X, INCX, WORK, INCX)
            ! WORK = alpha*work
            CALL ZSCAL(N, ALPHA, WORK, INCX)
            ! work = op(A)*work
            CALL ZTRMV(UPLOS(I),TRANS(J),DIAGS(K),N,A,N,WORK,INCX)
            ! work = beta*y + work
            CALL ZAXPY(N, BETA, Ys, 1, WORK, INCX)
            ! Compute work = -y + work
            CALL ZAXPY(N, NEG_ONE, Y, 1, WORK, INCX)
            ! Compute ||work||_2/||Y||_2
            NORM_ERR = DNRM2(N,WORK,INCX) / DNRM2(N,Y,1)
            WRITE(*,*) NORM_ERR
         END DO
         END DO
         END DO
         END DO
   10    DEALLOCATE(A)
         DEALLOCATE(As)
         DEALLOCATE(X)
         DEALLOCATE(Xs)
         DEALLOCATE(Y)
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
