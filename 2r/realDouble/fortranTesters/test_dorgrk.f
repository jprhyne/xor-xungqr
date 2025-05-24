      SUBROUTINE TEST_DORGRK(M, N)
*
*        Arguments
*
         INTEGER  M, N
*
*        Local variables
*
*        Scalars
*
         DOUBLE PRECISION  NORMA, NORM_ORTH, NORM_REPRES
         INTEGER           LWORK, I, J, INFO
*
*        Arrays
*
         DOUBLE PRECISION, ALLOCATABLE :: A(:,:), As(:,:), R(:,:),
     $      WORKMAT(:,:)

         DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TAU, WORK
*
*        Intrinsic Functions
*
         INTRINSIC         DBLE, SQRT
*
*        External Subroutines
*
         EXTERNAL          DGEMM, DGEQRF, DLACPY, DLARFT, DLASET,
     $                     DORGRK
*
*        External Functions
*
         DOUBLE PRECISION, EXTERNAL :: DLANGE
*
*        Parameters
*
         DOUBLE PRECISION ONE, ZERO, DNEG_ONE
         INTEGER           INEG_ONE
         PARAMETER(ONE=1.0D+0, ZERO=0.0D+0, INEG_ONE=-1.0D+0,
     $         DNEG_ONE=-1.0D+0)
*
*        Allocate our memory
*
         ALLOCATE(A(M,N))
         ALLOCATE(As(M,N))
         ALLOCATE(R(M,M))
         ALLOCATE(TAU(M))
         ALLOCATE(WORK(1))
         ALLOCATE(WORKMAT(M,M))
*
*        Set R to be the 0 matrix
*
         CALL DLASET('All', M, M, ZERO, ZERO, R, M)
*
*        Generate our data matrix A randomly
*
         CALL RANDOM_NUMBER(A)
*
*        Copy A into As for representation error checking
*
         CALL DLACPY('All', M, N, A, M, As, M)
*
*        Determine the size of work needed
*
         CALL DGERQF(M, N, A, M, TAU, WORK, INEG_ONE, INFO)
         LWORK = WORK(1)
         DEALLOCATE(WORK)
         ALLOCATE(WORK(LWORK))
*
*        Factorize A
*
         CALL DGERQF(M, N, A, M, TAU, WORK, LWORK, INFO)
*
*        Copy the upper triangular part of A (R) into R
*
         CALL DLACPY('Upper', M, M, A(1, N-M+1), M, R, M)
*
*        Compute the T matrix associated with V and Tau
*
         CALL MY_DLARFT_REC('Transpose', 'Rowwise', N, M, A, M, TAU,
     $            A(1, N-M+1), M)
*
*        Compute Q using our routine
*
         CALL DORGRK(M, N, A, M)
*
*        Determine if Q is orthogonal
*
*        First, we must set WORKMAT to be I_m
*
         CALL DLASET('All', M, M, ZERO, ONE, WORKMAT, M)
*
*        Next, we compute Q*Q' - I
*
         CALL DGEMM('No Transpose', 'Transpose', M, M, N, ONE, A,
     $         M, A, M, DNEG_ONE, WORKMAT, M)
*
*        Compute ||Q**T*Q - I_n||_F
*
         NORM_ORTH = ZERO
         DO I = 1, M
            DO J = 1, M
               NORM_ORTH = NORM_ORTH + WORKMAT(I,J)*WORKMAT(I,J)
            END DO
         END DO
*
*        Compute ||Q*Q' - I_n||_F / ||I_n||_F
*
         NORM_ORTH = SQRT(NORM_ORTH) / SQRT(DBLE(N))
*
*        Now, we compute the representation error
*
*        We need the norm of A for this
*
         NORMA = ZERO
         DO I = 1, M
            DO J = 1, N
               NORMA = NORMA + As(I,J)*As(I,J)
            END DO
         END DO
         NORMA = SQRT(NORMA)
*
*        Compute As = R*Q - As
*
         CALL DGEMM('No Transpose', 'No Transpose', M, N, M, ONE,
     $         R, M, A, M, DNEG_ONE, As, M)
*
*        Compute ||R*Q - As||_F
*
         NORM_REPRES = ZERO
         DO I = 1, M
            DO J = 1, N
               NORM_REPRES = NORM_REPRES + As(I,J)*As(I,J)
            END DO
         END DO
*
*        Compute ||R*Q - As||_F / ||As||_F
*
         NORM_REPRES = SQRT(NORM_REPRES) / NORMA
*
*        Print our error information
*
         WRITE(*,*) "representation norm: ", NORM_REPRES
         WRITE(*,*) "orthogonal norm:     ", NORM_ORTH
*
*        Free our memory
*
         DEALLOCATE(A)
         DEALLOCATE(As)
         DEALLOCATE(R)
         DEALLOCATE(TAU)
         DEALLOCATE(WORK)
         DEALLOCATE(WORKMAT)
      END SUBROUTINE
