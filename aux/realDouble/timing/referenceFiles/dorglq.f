      SUBROUTINE DORGLQ_REF(M,N,K,NB,A,LDA,TAU,WORK,LWORK,INFO)
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
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
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
      LWKOPT = MAX( 1, M )*NB
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
      NBMIN = 2
      NX = 0
      IWS = M
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'DORGLQ', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = M
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGLQ', ' ', M, N, K,
     $                      -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the last block.
*        The first kk rows are handled by the block method.
*
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
*
*        Set A(kk+1:m,1:kk) to zero.
*
         DO 20 J = 1, KK
            DO 10 I = KK + 1, M
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the last or only block.
*
      IF( KK.LT.M )
     $   CALL DORGL2( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,
     $                TAU( KK+1 ), WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.M ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL DLARFT( 'Forward', 'Rowwise', N-I+1, IB, A( I,
     $                      I ),
     $                      LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H**T to A(i+ib:m,i:n) from the right
*
               CALL DLARFB( 'Right', 'Transpose', 'Forward',
     $                      'Rowwise',
     $                      M-I-IB+1, N-I+1, IB, A( I, I ), LDA, WORK,
     $                      LDWORK, A( I+IB, I ), LDA, WORK( IB+1 ),
     $                      LDWORK )
            END IF
*
*           Apply H**T to columns i:n of current block
*
            CALL DORGL2( IB, N-I+1, IB, A( I, I ), LDA, TAU( I ),
     $                   WORK,
     $                   IINFO )
*
*           Set columns 1:i-1 of current block to zero
*
            DO 40 J = 1, I - 1
               DO 30 L = I, I + IB - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGLQ_REF
*
      END SUBROUTINE
      !
      SUBROUTINE DORGLQ_LVL2(M,N,K,NB,A,LDA,TAU, WORK, LWORK, INFO )
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
      NBMIN = 2
      NX = MAX( 0, ILAENV( 3, 'DORGLQ', ' ', M, N, K, -1 ) )
      IWS = M
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Handle the first block assuming we are applying to the
*        identity, then resume regular blocking method after
*
         KI = K - 2 * NB
         KK = K - NB
      ELSE
         KK = 0
      END IF
*
*     Potentially bail to the unblocked version
*
      IF( KK.EQ.0 ) THEN
         CALL DORGL2( M, N, K, A, LDA, TAU, WORK, IINFO )
      END IF
*
      IF( KK.GT.0 ) THEN
         I = KK + 1
         IB = NB
*
*        Form the triangular factor of the block reflector
*        H = H(i) H(i+1) . . . H(i+ib-1)
*
         CALL DLARFT_LVL2( 'Forward', 'Transpose', N-I+1, IB, A( I, I ),
     $               LDA, TAU( I ), A( I, I ), LDA )
*
*        Apply H to A(i+ib:m,i:n) from the right
*
         CALL DLARFB0C2('Identity', 'B', 'Right', 'No Transpose',
     $         'Forward', 'Rowwise', M-I-IB+1, N-I+1, IB, A(I,I), LDA,
     $         A(I,I), LDA, A(I+IB,I), LDA)
*
*        Apply H to columns i:n of current block

         CALL DORGLK( IB, N-I+1, A( I, I ), LDA)
*
*        Use blocked code
*
         DO I = KI + 1, 1, -NB
            IB = NB
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL DLARFT_LVL2( 'Forward', 'Transpose', N-I+1, IB, A(I,I),
     $                  LDA, TAU( I ), A( I, I ), LDA )
*
*           Apply H to A(i+ib:m,i:n) from the right
*
            CALL DLARFB0C2('General', 'B', 'Right', 'No Transpose',
     $            'Forward', 'Rowwise', M-I-IB+1, N-I+1, IB, A(I,I),
     $            LDA, A(I,I), LDA, A(I+IB,I), LDA)
*
*           Apply H to columns i:n of current block
*
            CALL DORGLK( IB, N-I+1, A( I, I ), LDA)
         END DO
*
*        This checks for if K was a perfect multiple of NB
*        so that we only have a special case for the last block when
*        necessary
*
         IF(I.LT.1) THEN
            IB = I + NB - 1
            I = 1
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL DLARFT_LVL2( 'Forward', 'Transpose', N-I+1, IB, A(I,I),
     $                  LDA, TAU( I ), A( I, I ), LDA )
*
*           Apply H to A(i+ib:m,i:n) from the right
*
            CALL DLARFB0C2('General', 'B', 'Right', 'No Transpose',
     $            'Forward', 'Rowwise', M-I-IB+1, N-I+1, IB, A(I,I),
     $            LDA, A(I,I), LDA, A(I+IB,I), LDA)
*
*           Apply H to columns i:n of current block
*
            CALL DORGLK( IB, N-I+1, A( I, I ), LDA)
         END IF
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGLQ_LVL2
*
      END SUBROUTINE
      !
      SUBROUTINE DORGLQ_REC(M,N,K,NB, A, LDA, TAU, WORK, LWORK, INFO)
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
      NBMIN = 2
      NX = MAX( 0, ILAENV( 3, 'DORGLQ', ' ', M, N, K, -1 ) )
      IWS = M
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Handle the first block assuming we are applying to the
*        identity, then resume regular blocking method after
*
         KI = K - 2 * NB
         KK = K - NB
      ELSE
         KK = 0
      END IF
*
*     Potentially bail to the unblocked version
*
      IF( KK.EQ.0 ) THEN
         CALL DORGL2( M, N, K, A, LDA, TAU, WORK, IINFO )
      END IF
*
      IF( KK.GT.0 ) THEN
         I = KK + 1
         IB = NB
*
*        Form the triangular factor of the block reflector
*        H = H(i) H(i+1) . . . H(i+ib-1)
*
         CALL DLARFT_REC( 'Forward', 'Transpose', N-I+1, IB, A( I, I ),
     $               LDA, TAU( I ), A( I, I ), LDA )
*
*        Apply H to A(i+ib:m,i:n) from the right
*
         CALL DLARFB0C2('Identity', 'B', 'Right', 'No Transpose',
     $         'Forward', 'Rowwise', M-I-IB+1, N-I+1, IB, A(I,I), LDA,
     $         A(I,I), LDA, A(I+IB,I), LDA)
*
*        Apply H to columns i:n of current block

         CALL DORGLK( IB, N-I+1, A( I, I ), LDA)
*
*        Use blocked code
*
         DO I = KI + 1, 1, -NB
            IB = NB
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL DLARFT_REC( 'Forward', 'Transpose', N-I+1, IB, A(I,I),
     $                  LDA, TAU( I ), A( I, I ), LDA )
*
*           Apply H to A(i+ib:m,i:n) from the right
*
            CALL DLARFB0C2('General', 'B', 'Right', 'No Transpose',
     $            'Forward', 'Rowwise', M-I-IB+1, N-I+1, IB, A(I,I),
     $            LDA, A(I,I), LDA, A(I+IB,I), LDA)
*
*           Apply H to columns i:n of current block
*
            CALL DORGLK( IB, N-I+1, A( I, I ), LDA)
         END DO
*
*        This checks for if K was a perfect multiple of NB
*        so that we only have a special case for the last block when
*        necessary
*
         IF(I.LT.1) THEN
            IB = I + NB - 1
            I = 1
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL DLARFT_REC( 'Forward', 'Transpose', N-I+1, IB, A(I,I),
     $                  LDA, TAU( I ), A( I, I ), LDA )
*
*           Apply H to A(i+ib:m,i:n) from the right
*
            CALL DLARFB0C2('General', 'B', 'Right', 'No Transpose',
     $            'Forward', 'Rowwise', M-I-IB+1, N-I+1, IB, A(I,I),
     $            LDA, A(I,I), LDA, A(I+IB,I), LDA)
*
*           Apply H to columns i:n of current block
*
            CALL DORGLK( IB, N-I+1, A( I, I ), LDA)
         END IF
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGLQ_REC
*
      END SUBROUTINE
      !
      SUBROUTINE DORGLQ_UT_INV(M,N,K,NB, A, LDA, TAU, WORK, LWORK, INFO)
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
      NBMIN = 2
      NX = MAX( 0, ILAENV( 3, 'DORGLQ', ' ', M, N, K, -1 ) )
      IWS = M
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Handle the first block assuming we are applying to the
*        identity, then resume regular blocking method after
*
         KI = K - 2 * NB
         KK = K - NB
      ELSE
         KK = 0
      END IF
*
*     Potentially bail to the unblocked version
*
      IF( KK.EQ.0 ) THEN
         CALL DORGL2( M, N, K, A, LDA, TAU, WORK, IINFO )
      END IF
*
      IF( KK.GT.0 ) THEN
         I = KK + 1
         IB = NB
*
*        Form the triangular factor of the block reflector
*        H = H(i) H(i+1) . . . H(i+ib-1)
*
         CALL DLARFT_UT( 'F','T', '1', N-I+1, IB, A( I, I ),
     $               LDA, TAU( I ), A( I, I ), LDA )
*
*        Apply H to A(i+ib:m,i:n) from the right
*
         CALL DLARFB0C2('Identity', 'B', 'Right', 'No Transpose',
     $         'Forward', 'Rowwise', M-I-IB+1, N-I+1, IB, A(I,I), LDA,
     $         A(I,I), LDA, A(I+IB,I), LDA)
*
*        Apply H to columns i:n of current block

         CALL DORGLK( IB, N-I+1, A( I, I ), LDA)
*
*        Use blocked code
*
         DO I = KI + 1, 1, -NB
            IB = NB
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL DLARFT_UT( 'F', 'T', '1', N-I+1, IB, A(I,I),
     $                  LDA, TAU( I ), A( I, I ), LDA )
*
*           Apply H to A(i+ib:m,i:n) from the right
*
            CALL DLARFB0C2('General', 'B', 'Right', 'No Transpose',
     $            'Forward', 'Rowwise', M-I-IB+1, N-I+1, IB, A(I,I),
     $            LDA, A(I,I), LDA, A(I+IB,I), LDA)
*
*           Apply H to columns i:n of current block
*
            CALL DORGLK( IB, N-I+1, A( I, I ), LDA)
         END DO
*
*        This checks for if K was a perfect multiple of NB
*        so that we only have a special case for the last block when
*        necessary
*
         IF(I.LT.1) THEN
            IB = I + NB - 1
            I = 1
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL DLARFT_UT( 'F', 'T', '1', N-I+1, IB, A(I,I),
     $                  LDA, TAU( I ), A( I, I ), LDA )
*
*           Apply H to A(i+ib:m,i:n) from the right
*
            CALL DLARFB0C2('General', 'B', 'Right', 'No Transpose',
     $            'Forward', 'Rowwise', M-I-IB+1, N-I+1, IB, A(I,I),
     $            LDA, A(I,I), LDA, A(I+IB,I), LDA)
*
*           Apply H to columns i:n of current block
*
            CALL DORGLK( IB, N-I+1, A( I, I ), LDA)
         END IF
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGLQ_UT_INV
*
      END SUBROUTINE
      !
      SUBROUTINE DORGLQ_UT_SOLVE(M,N,K,NB,A,LDA,TAU,WORK,LWORK,INFO)
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
      NBMIN = 2
      NX = MAX( 0, ILAENV( 3, 'DORGLQ', ' ', M, N, K, -1 ) )
      IWS = M
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Handle the first block assuming we are applying to the
*        identity, then resume regular blocking method after
*
         KI = K - 2 * NB
         KK = K - NB
      ELSE
         KK = 0
      END IF
*
*     Potentially bail to the unblocked version
*
      IF( KK.EQ.0 ) THEN
         CALL DORGL2( M, N, K, A, LDA, TAU, WORK, IINFO )
      END IF
*
      IF( KK.GT.0 ) THEN
         I = KK + 1
         IB = NB
*
*        Form the triangular factor of the block reflector
*        H = H(i) H(i+1) . . . H(i+ib-1)
*
         CALL DLARFT_UT( 'F','T', '2', N-I+1, IB, A( I, I ),
     $               LDA, TAU( I ), A( I, I ), LDA )
*
*        Apply H to A(i+ib:m,i:n) from the right
*
         CALL DLARFB0C2('Identity', 'A', 'Right', 'No Transpose',
     $         'Forward', 'Rowwise', M-I-IB+1, N-I+1, IB, A(I,I), LDA,
     $         A(I,I), LDA, A(I+IB,I), LDA)
*
*        Apply H to columns i:n of current block

         CALL DORGLK( IB, N-I+1, A( I, I ), LDA)
*
*        Use blocked code
*
         DO I = KI + 1, 1, -NB
            IB = NB
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL DLARFT_UT( 'F', 'T', '2', N-I+1, IB, A(I,I),
     $                  LDA, TAU( I ), A( I, I ), LDA )
*
*           Apply H to A(i+ib:m,i:n) from the right
*
            CALL DLARFB0C2('General', 'A', 'Right', 'No Transpose',
     $            'Forward', 'Rowwise', M-I-IB+1, N-I+1, IB, A(I,I),
     $            LDA, A(I,I), LDA, A(I+IB,I), LDA)
*
*           Apply H to columns i:n of current block
*
            CALL DORGLK( IB, N-I+1, A( I, I ), LDA)
         END DO
*
*        This checks for if K was a perfect multiple of NB
*        so that we only have a special case for the last block when
*        necessary
*
         IF(I.LT.1) THEN
            IB = I + NB - 1
            I = 1
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL DLARFT_UT( 'F', 'T', '2', N-I+1, IB, A(I,I),
     $                  LDA, TAU( I ), A( I, I ), LDA )
*
*           Apply H to A(i+ib:m,i:n) from the right
*
            CALL DLARFB0C2('General', 'A', 'Right', 'No Transpose',
     $            'Forward', 'Rowwise', M-I-IB+1, N-I+1, IB, A(I,I),
     $            LDA, A(I,I), LDA, A(I+IB,I), LDA)
*
*           Apply H to columns i:n of current block
*
            CALL DORGLK( IB, N-I+1, A( I, I ), LDA)
         END IF
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGLQ_UT_INV
*
      END SUBROUTINE
