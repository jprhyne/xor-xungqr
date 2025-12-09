      SUBROUTINE DORGLQ_REF_K(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,
     $                   LWKOPT, NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORGL2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LWKOPT = 2*MAX( 1, M )*K
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGLQ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
      LDWORK = 2*M
*
*     Set A(k+1:m,1:k) to zero.
*
      CALL DLASET('All', M-K, K, ZERO, ZERO, A(K+1,1), LDA)
*
*     Set A(k+1:m,k+1:n) to eye(m-k,n-k)
*
      CALL DLASET('All', M-K, N-K, ZERO, ONE, A(K+1,K+1), LDA)
*
*     Form the triangular factor of the block reflector
*     H = H(1) H(2) . . . H(k)
*
      CALL DLARFT('Forward','Rowwise',N,K,A,LDA,TAU,WORK,LDWORK)
*
*     Apply H**T to A(k+1:m,1:n) from the right
*
      CALL DLARFB('Right','Transpose','Forward','Rowwise',M-K,N,K,
     $             A,LDA,WORK,LDWORK,A(K+1,1),LDA,WORK(K+1),LDWORK)
*
*     Apply H**T to A(1:k,1:n)
*
      CALL DORGL2( K,N,K,A,LDA,TAU,WORK,IINFO )
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGLQ_REF
*
      END SUBROUTINE
      !
      SUBROUTINE DORGLQ_LVL2_K(M,N,K,A,LDA,TAU, WORK, LWORK, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, KI, KK, LWKOPT,
     $                   NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB0C2, DLARFT_LVL2, DORGL2,
     $                   DORGLK, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      ! Only need a workspace for calls to dorgl2
      LWKOPT = MAX( 1, M )
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGLQ_LVL2', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*
*     Form the triangular factor of the block reflector
*     H = H(1) H(2) . . . H(k)
*
      CALL DLARFT_LVL2('Forward','Transpose',N,K,A,LDA,TAU,A,LDA)
*
*     Apply H to A(k+1:m,1:n) from the right
*
      CALL DLARFB0C2('Identity', 'B', 'Right', 'No Transpose',
     $      'Forward', 'Rowwise', M-K, N, K, A, LDA, A, LDA,
     $      A(K+1,1), LDA)
*
*     Apply H to A(1:k,1:n)

      CALL DORGLK(K,N,A,LDA)
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGLQ_LVL2
*
      END SUBROUTINE
      !
      SUBROUTINE DORGLQ_REC_K(M,N,K, A, LDA, TAU, WORK, LWORK, INFO)
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, KI, KK, LWKOPT,
     $                   NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB0C2, DLARFT_REC, DORGL2,
     $                   DORGLK, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      ! Only need a workspace for calls to dorgl2
      LWKOPT = MAX( 1, M )
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGLQ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*
*     Form the triangular factor of the block reflector
*     H = H(1) H(2) . . . H(k)
*
      CALL DLARFT_REC('Forward','Transpose',N,K,A,LDA,TAU,A,LDA)
*
*     Apply H to A(k+1:m,1:n) from the right
*
      CALL DLARFB0C2('Identity', 'B', 'Right', 'No Transpose',
     $      'Forward', 'Rowwise', M-K, N, K, A, LDA, A, LDA,
     $      A(K+1,1), LDA)
*
*     Apply H to A(1:k,1:n)

      CALL DORGLK(K,N,A,LDA)
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGLQ_REC
*
      END SUBROUTINE
      !
      SUBROUTINE DORGLQ_UT_INV_K(M,N,K, A, LDA, TAU, WORK, LWORK, INFO)
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, KI, KK, LWKOPT,
     $                   NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB0C2, DLARFT_UT, DORGL2,
     $                   DORGLK, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      ! Only need a workspace for calls to dorgl2
      LWKOPT = MAX( 1, M )
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGLQ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*
*     Form the triangular factor of the block reflector
*     H = H(1) H(2) . . . H(k)
*
      CALL DLARFT_UT('Forward','Transpose','1',N,K,A,LDA,TAU,A,LDA)
*
*     Apply H to A(k+1:m,1:n) from the right
*
      CALL DLARFB0C2('Identity', 'B', 'Right', 'No Transpose',
     $      'Forward', 'Rowwise', M-K, N, K, A, LDA, A, LDA,
     $      A(K+1,1), LDA)
*
*     Apply H to A(1:k,1:n)

      CALL DORGLK(K,N,A,LDA)
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGLQ_UT_INV
*
      END SUBROUTINE
      !
      SUBROUTINE DORGLQ_UT_SOLVE_K(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, KI, KK, LWKOPT,
     $                   NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB0C2, DLARFT_UT, DORGL2,
     $                   DORGLK, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      ! Only need a workspace for calls to dorgl2
      LWKOPT = MAX( 1, M )
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGLQ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*
*     Form the triangular factor of the block reflector
*     H = H(1) H(2) . . . H(k)
*
      CALL DLARFT_UT('Forward','Transpose','2',N,K,A,LDA,TAU,A,LDA)
*
*     Apply H to A(k+1:m,1:n) from the right
*
      CALL DLARFB0C2('Identity', 'A', 'Right', 'No Transpose',
     $      'Forward', 'Rowwise', M-K, N, K, A, LDA, A, LDA,
     $      A(K+1,1), LDA)
*
*     Apply H to A(1:k,1:n)

      CALL DORGLK(K,N,A,LDA)
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGLQ_UT_INV
*
      END SUBROUTINE
