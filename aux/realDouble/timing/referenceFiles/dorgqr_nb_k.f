      SUBROUTINE DORGQR_REF_K(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
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
     $                   LWKOPT, NB, NBMIN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORG2R, XERBLA
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
      LWKOPT = 2*MAX( 1, N )*K
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.LWKOPT .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
      LDWORK = 2*N
*
*     Set A(1:k, k+1:n) to zero
*
      CALL DLASET('All', K, N-K, ZERO, ZERO, A(1,K+1), LDA)
*
*     Set A(k+1:m, k+1:n) to eye(m-k,n-k)
*
      CALL DLASET('All', M-K, N-K, ZERO, ONE, A(K+1,K+1), LDA)
*
*     Form the triangular factor of the block reflector
*     H = H(1) H(2) . . . H(k)
*
      CALL DLARFT('Forward', 'Columnwise', M, K, A, LDA, TAU,
     $            WORK, LDWORK)
*
*     Apply H to A(1:m,k+1:n) from the left
*
      CALL DLARFB( 'Left', 'No transpose', 'Forward',
     $             'Columnwise', M, N-K, K, A, LDA, WORK, LDWORK,
     $             A( 1, K+1 ), LDA, WORK( K+1 ), LDWORK )
*
*     Apply H to rows i:m of current block
*
      CALL DORG2R( M, K, K, A, LDA, TAU, WORK, IINFO )
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGQR
*
      END SUBROUTINE
*
      SUBROUTINE DORGQR_OPT_K(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
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
      INTEGER            I, IB, IINFO, KI, KK, LWKOPT,
     $                   NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL             DLARFB0C2, DLARFT, DORG2R,
     $                     DORGKR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
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
*
*     Only need a workspace for dorg2r in case of bail out
*
      LWKOPT = MAX( 1, N )
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQR_REC', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*
*     Form the triangular factor of the block reflector
*     H = H(1) H(2) . . . H(i+ib-1)
*
      CALL DLARFT('Forward', 'Column', M, K, A, LDA, TAU, A, LDA)
*
*     Apply H to A(i:m,i+ib:n) from the left
*
      CALL DLARFB0C2('Identity', 'B', 'Left', 'No Transpose',
     $   'Forward', 'Column', M, N-K, K, A, LDA, A, LDA, A(1,K+1), LDA)
*
*     Apply H to rows i:m of current block
*
      CALL DORGKR('1',M, K, A, LDA)
*
      WORK( 1 ) = N
      RETURN
*
*     End of DORGQR_REC
*
      END SUBROUTINE
      SUBROUTINE DORGQR_REC_K(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
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
      INTEGER            I, IB, IINFO, KI, KK, LWKOPT,
     $                   NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL             DLARFB0C2, DLARFT_REC, DORG2R,
     $                     DORGKR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
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
*
*     Only need a workspace for dorg2r in case of bail out
*
      LWKOPT = MAX( 1, N )
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQR_REC', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*
*     Form the triangular factor of the block reflector
*     H = H(1) H(2) . . . H(i+ib-1)
*
      CALL DLARFT_REC('Forward', 'Column', M, K, A, LDA, TAU, A, LDA)
*
*     Apply H to A(i:m,i+ib:n) from the left
*
      CALL DLARFB0C2('Identity', 'B', 'Left', 'No Transpose',
     $   'Forward', 'Column', M, N-K, K, A, LDA, A, LDA, A(1,K+1), LDA)
*
*     Apply H to rows i:m of current block
*
      CALL DORGKR('1',M, K, A, LDA)
*
      WORK( 1 ) = N
      RETURN
*
*     End of DORGQR_REC
*
      END SUBROUTINE

      SUBROUTINE DORGQR_LVL2_K(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
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
      INTEGER            I, IB, IINFO, KI, KK, LWKOPT,
     $                   NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL             DLARFB0C2, DLARFT_LVL2, DORG2R,
     $                     DORGKR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
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
*
*     Only need a workspace for dorg2r in case of bail out
*
      LWKOPT = MAX( 1, N )
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQR_LVL2', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*
*     Form the triangular factor of the block reflector
*     H = H(1) H(2) . . . H(i+ib-1)
*
      CALL DLARFT_LVL2('Forward', 'Column', M, K, A, LDA, TAU, A, LDA)
*
*     Apply H to A(i:m,i+ib:n) from the left
*
      CALL DLARFB0C2('Identity', 'B', 'Left', 'No Transpose',
     $   'Forward', 'Column', M, N-K, K, A, LDA, A, LDA, A(1,K+1), LDA)
*
*     Apply H to rows i:m of current block
*
      CALL DORGKR('1',M, K, A, LDA)
*
      WORK( 1 ) = N
      RETURN
*
*     End of DORGQR_LVL2
*
      END SUBROUTINE
*
      SUBROUTINE DORGQR_UT_INV_K(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
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
      INTEGER            I, IB, IINFO, KI, KK, LWKOPT,
     $                   NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL             DLARFB0C2, DLARFT_UT, DORG2R,
     $                     DORGKR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
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
*
*     Only need a workspace for dorg2r in case of bail out
*
      LWKOPT = MAX( 1, N )
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQR_UT_INV', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*
*     Form the triangular factor of the block reflector
*     H = H(1) H(2) . . . H(i+ib-1)
*
      CALL DLARFT_UT('Forward','Column','1',M,K,A,LDA,TAU,A,LDA)
*
*     Apply H to A(i:m,i+ib:n) from the left
*
      CALL DLARFB0C2('Identity', 'B', 'Left', 'No Transpose',
     $   'Forward', 'Column', M, N-K, K, A, LDA, A, LDA, A(1,K+1), LDA)
*
*     Apply H to rows i:m of current block
*
      CALL DORGKR('1',M, K, A, LDA)
*
      WORK( 1 ) = N
      RETURN
*
*     End of DORGQR_UT_INV
*
      END SUBROUTINE
      SUBROUTINE DORGQR_UT_SOLVE_K(M,N,K,A,LDA,TAU,WORK,LWORK,INFO)
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
      INTEGER            I, IB, IINFO, KI, KK, LWKOPT,
     $                   NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL             DLARFB0C2, DLARFT_UT, DORG2R,
     $                     DORGKR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
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
*
*     Only need a workspace for dorg2r in case of bail out
*
      LWKOPT = MAX( 1, N )
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQR_UT_SOLVE', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*
*     Form the triangular factor of the block reflector
*     H = H(1) H(2) . . . H(i+ib-1)
*
      CALL DLARFT_UT('Forward','Column','2',M,K,A,LDA,TAU,A,LDA)
*
*     Apply H to A(i:m,i+ib:n) from the left
*
      CALL DLARFB0C2('Identity', 'A', 'Left', 'No Transpose',
     $   'Forward', 'Column', M, N-K, K, A, LDA, A, LDA, A(1,K+1), LDA)
*
*     Apply H to rows i:m of current block
*
      CALL DORGKR('1',M, K, A, LDA)
*
      WORK( 1 ) = N
      RETURN
*
*     End of DORGQR_UT_SOLVE
*
      END SUBROUTINE
