      PROGRAM TEST_DLUSM
      INTEGER              N, I, J, K, L
      DOUBLE PRECISION     NORM_ERR, NORM_X, ONE, ZERO, ALPHA,
      CHARACTER            SIDE, UPLOA, DIAGA, DIAGB
      DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:), X(:,:), WORK(:),
     $               WORKMAT(:,:)
      CHARACTER :: SIDES(2), UPLOAS(2), DIAGAS(2), DIAGBS(2)
      DOUBLE PRECISION  DLANGE
      EXTERNAL          DLANGE

      READ(*,*) N
      ONE = 1.0D+0
      ZERO = 0.0D+0

      ALLOCATE(A(N,N))
      ALLOCATE(B(N,N))
      ALLOCATE(X(N,N))
      ALLOCATE(WORKMAT(N,N))
      ALLOCATE(WORK(N))

      DO I = 1, 2 ! SIDES = ('L', 'R')
      DO J = 1, 2 ! UPLOAS = ('U', 'L')
      DO K = 1, 2 ! DIAGAS = ('U', 'N')
      DO L = 1, 2 ! DIAGBS = ('U', 'N')
         ! Generate our matrices
      END DO
      END DO
      END DO
      END DO
      END PROGRAM
