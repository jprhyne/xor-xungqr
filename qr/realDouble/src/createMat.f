*     How is this different from just using DLAGGE?
*     Optionally allow a seed?
      SUBROUTINE CREATE_DMAT( M, N, SIGMA, A, LDA, INFO )
*     Input arguments
      INTEGER           M, N, LDA, INFO, LWORK
      DOUBLE PRECISION  SIGMA( M ), A( LDA, *)

*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )

*
*     Local Arrays
*
      DOUBLE PRECISION  U(M, M), V(N, N), S(M, N), TAU1(M), TAU2(N)
      DOUBLE PRECISION  W(M,N), WORKARR(2)
*     allocatable array to be used for workspaces
      real*8, allocatable :: WORK(:)
*
*     Local Scalars
*
      INTEGER           I, J, NB
*
*     External subroutines
*
      EXTERNAL          MY_DORGQR, DEQGRF
*
*     Instrinsic Functions
*
      INTRINSIC         RAND, MAX
*
*     Testing input arguments
*
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M) THEN
         INFO = -2
      ELSE IF( LDA.GT.MAX(1,M) ) THEN
         INFO = -3
      END IF

      IF( INFO.NE.0 ) THEN
         RETURN
      END IF

*     Create S = DIAG(S)
      DO 20 J = 1, N
         DO 10 I = 1, J-1
            S(I,J) = 0
   10    CONTINUE
         S(J,J) = SIGMA(J)
         DO 15 I = J+1, M
            S(I,J) = 0
   15    CONTINUE
   20 CONTINUE
*
*     Create U \in \R^{m by m} orthogonal
*

*
*     Initialize U as a random matrix
*
      DO 40 J = 1, M
         DO 30 I = 1,M
            U(I,J) = RAND()
   30    CONTINUE
   40 CONTINUE

*
*     Query the workspace
*
      LWORK=-1
      CALL DGEQRF(M,M,U,M,TAU1,WORKARR,LWORK,INFO)
      CALL MY_DORGQR(M,M,M,NB,U,M,TAU1,WORKARR(2),LWORK, INFO)
*
*     Allocate the workspace
*
      LWORK = MAX(WORKARR(1),WORKARR(2))
      ALLOCATE(WORK(LWORK))
*
*     Construct U to be orthogonal
*
      CALL DGEQRF(M,M,U,M,TAU1,WORK,LWORK,INFO)
      CALL MY_DORGQR(M,M,M,NB,U,M,TAU1,WORK,LWORK,INFO)

*
*     Create V \in \R^{n by n} orthogonal
*
      DO 60 J = 1, N
         DO 50 I = 1,N
            V(I,J) = RAND()
   50    CONTINUE
   60 CONTINUE
*
*     Query the workspace
*
      LWORK = -1
      CALL DGEQRF(N,N,V,N,TAU2,WORKARR,LWORK,INFO)
      CALL MY_DORGQR(N,N,N,NB,V,M,TAU2,WORKARR(2),LWORK, INFO)
*
*     Allocate the workspace
*
      LWORK = MAX(WORKARR(1),WORKARR(2))
      DEALLOCATE(WORK)
      ALLOCATE(WORK(LWORK))
*
*     Construct V to be orthogonal
*
      CALL DGEQRF(N,N,V,N,TAU2,WORK,LWORK,INFO)
      CALL MY_DORGQR(N,N,N,NB,V,N,TAU2,WORK,LWORK,INFO)
      DEALLOCATE(WORK)
*
* Compute W = U * \Sigma
*
      CALL DGEMM('N', 'N', M, M, N, ONE, U, M, S, M, ZERO, W, M)
*
* Compute A = W * V**T = U * \Sigma * V**T
*
      CALL DGEMM('N', 'T', M, N, N, ONE, W, M, V, N, ZERO, A, LDA)
      END
