*> \brief \b DTRTRM computes an in place triangular-triangular matrix product
*
*  =========== DOCUMENTATION ===========
*
*  Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     RECURSIVE SUBROUTINE DTRTRM(SIDE, UPLO, TRANSV, DIAGT, DIAGV,
*    $                        N, ALPHA, T, LDT, V, LDV)
*
*        .. Scalar Arguments ..
*        INTEGER           N, LDT, LDV
*        CHARACTER         SIDE, UPLO, TRANSV, DIAGT, DIAGV
*        DOUBLE PRECISION  ALPHA
*        ..
*        .. Array Arguments ..
*        DOUBLE PRECISION  T(LDT,*), V(LDV,*)
*        ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DTRMMOOP performs one  of the matrix-matrix operatiions
*>
*>       T = \alpha op(V) * T
*>                      or
*>       T = \alpha T * op(V)
*> where \alpha is a scalar, T and V are unit, or non-unit, upper or
*> lower triangular matrix, and op(V) is one of
*>
*>       op(V) = V      or       op(V) = V**T
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>           On entry, SIDE specifies whether op(V) multiplies T from
*>           the left or right as follows:
*>
*>             SIDE = 'L' or 'l'    T = \alpha op(V) * T
*>
*>             SIDE = 'R' or 'r'    T = \alpha T * op(V)
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix T is an upper or
*>           lower triangular matrix as follows:
*>             UPLO = 'U' or 'u'    T is upper triangular
*>
*>             UPLO = 'L' or 'l'    T is lower triangular
*> \Endverbatim
*>
*> \param[in] TRANSV
*> \verbatim
*>          TRANSV is CHARACTER*1
*>           On entry, TRANSV specifies the form of op(V) to be used in
*>           the matrix multiplication as follows:
*>             TRANSV = 'N' or 'n'    op(V) = V
*>
*>             TRANSV = 'T' or 't'    op(V) = V**T
*>
*>             TRANSV = 'C' or 'c'    op(V) = V**T
*> \endverbatim
*>
*> \param[in] DIAGT
*> \verbatim
*>          DIAGT is CHARACTER*1
*>           On entry, DIAGT specifies whether or not T is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u'      T is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'      T is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] DIAGV
*> \verbatim
*>          DIAGV is CHARACTER*1
*>           On entry, DIAGV specifies whether or not V is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u'      V is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'      V is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of rows and columns of T.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha. When alpha is
*>           zero then T and V are not referenced, and T and V need not
*>           be set before entry.
*> \endverbatim
*>
*> \param[in] T
*> \verbatim
*>          T is DOUBLE PRECISION array, dimension ( LDT, N )
*>           Before entry with UPLO = 'U' or 'u', the leading k-by-k
*>           upper triangular part of the array T must contain the upper
*>           triangular matrix and the strictly lower triangular part of
*>           T is not referenced.
*>           Before entry  with  UPLO = 'L' or 'l', the leading k-by-k
*>           lower triangular part of the array T must contain the lower
*>           triangular matrix and the strictly upper triangular part of
*>           T is not referenced.
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*>           T  are not referenced either,  but are assumed to be  unity.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>           On entry, LDT specifies the first dimension of T as declared
*>           in the calling (sub) program. LDT must be at least max( 1, n ).
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is DOUBLE PRECISION array, dimension ( LDV, N )
*>           Before entry with UPLO = 'U' or 'u', the leading k-by-k
*>           upper triangular part of the array V must contain the upper
*>           triangular matrix and the strictly lower triangular part of
*>           V is not referenced.
*>           Before entry  with  UPLO = 'L' or 'l', the leading k-by-k
*>           lower triangular part of the array V must contain the lower
*>           triangular matrix and the strictly upper triangular part of
*>           V is not referenced.
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*>           V  are not referenced either,  but are assumed to be  unity.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>           On entry, LDV specifies the first dimension of T as declared
*>           in the calling (sub) program. LDV must be at least max( 1, n ).
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION.
*>           On entry, BETA specifies the scalar beta. When beta is
*>           zero then C is not referenced on entry, and C need not
*>           be set before entry.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension ( LDB, N )
*>           Before entry, the leading m-by-n part of the array C must
*>           contain the matrix C, and on exit is overwritten by the
*>           transformed matrix.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>           On entry, LDC specifies the first dimension of C as declared
*>           in the calling (sub) program. LDC must be at least
*>           max( 1, m ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*  =====================================================================
      RECURSIVE SUBROUTINE DTRTRM(SIDE, UPLO, TRANSV, DIAGT, DIAGV,
     $                        N, ALPHA, T, LDT, V, LDV)
*
*        .. Scalar Arguments ..
         INTEGER           N, LDT, LDV
         CHARACTER         SIDE, UPLO, TRANSV, DIAGT, DIAGV
         DOUBLE PRECISION  ALPHA
*        ..
*        .. Array Arguments ..
         DOUBLE PRECISION  T(LDT,*), V(LDV,*)
*        ..
*
*  =====================================================================
*
*        .. External Functions ..
         LOGICAL           LSAME
         EXTERNAL          LSAME
*        ..
*        .. External Subroutines ..
         EXTERNAL          DTRMM, DTRMMOOP
*        ..
*        .. Local Scalars ..
         INTEGER           K, INFO
         LOGICAL           TLEFT, TUPPER, VTRANS, VUNIT, TUNIT
*        ..
*        .. Local Parameters ..
         DOUBLE PRECISION ONE
         PARAMETER(ONE=1.0D+0)
*        ..
*
*        Beginning of Executable Statements
*
         TUNIT = LSAME(DIAGT, 'U')
         VUNIT = LSAME(DIAGV, 'U')
*
*        Terminating Case
*
         IF (N.EQ.1) THEN
         ELSE IF(N.LE.0) THEN
            RETURN
         END IF
*
*        Recursive case
*
         TUPPER = LSAME(UPLO,'U')
         TLEFT  = LSAME(SIDE,'L')
         VTRANS = LSAME(TRANSV,'T').OR.LSAME(TRANSV,'C')

         K = N / 2
         IF(TUPPER) THEN
            IF(TLEFT) THEN
               IF(VTRANS) THEN
                  CALL DTRMM('Right', 'Lower', TRANSV, DIAGV, K,
     $                     N-K, ALPHA, V(K+1, K+1), LDV, T(1, K+1), LDT)
                  CALL DTRMMOOP(SIDE, UPLO, 'No Transpose', TRANSV,
     $                     DIAGT, K, N-K, ALPHA, T, LDT, V(K+1, 1), LDV,
     $                     ONE, T(1, K+1), LDT)
               ELSE
                  CALL DTRMM('Right', 'Upper', TRANSV, DIAGV, K,
     $                     N-K, ALPHA, V(K+1, K+1), LDV, T(1, K+1), LDT)
                  CALL DTRMMOOP(SIDE, UPLO, 'No Transpose', TRANSV,
     $                     DIAGT, K, N-K, ALPHA, T, LDT, V(1, K+1), LDV,
     $                     ONE, T(1, K+1), LDT)
               END IF
            ELSE
               IF(VTRANS) THEN
                  CALL DTRMM('Left', 'Lower', TRANSV, DIAGV, K,
     $                     N-K, ALPHA, V, LDV, T(1, K+1), LDT)
                  CALL DTRMMOOP(SIDE, UPLO, 'No Transpose', TRANSV,
     $                     DIAGT, K, N-K, ALPHA, T(K+1, K+1), LDT,
     $                     V(K+1, 1), LDV, ONE, T(1, K+1), LDT)
               ELSE
                  CALL DTRMM('Left', 'Upper', TRANSV, DIAGV, K,
     $                     N-K, ALPHA, V, LDV, T(1, K+1), LDT)
                  CALL DTRMMOOP(SIDE, UPLO, 'No Transpose', TRANSV,
     $                     DIAGT, K, N-K, ALPHA, T(K+1, K+1), LDT,
     $                     V(1, K+1), LDV, ONE, T(1, K+1), LDT)
               END IF
            END IF
         ELSE
            IF(TLEFT) THEN
               IF(VTRANS) THEN
                  CALL DTRMM('Right', 'Upper', TRANSV, DIAGV, N-K,
     $                     K, ALPHA, V, LDV, T(K+1, 1), LDT)
                  CALL DTRMMOOP(SIDE, UPLO, 'No Transpose', TRANSV,
     $                     DIAGT, N-K, K, ALPHA, T(K+1, K+1), LDT,
     $                     V(1, K+1), LDV, ONE, T(K+1, 1), LDT)
               ELSE
                  CALL DTRMM('Right', 'Lower', TRANSV, DIAGV, N-K,
     $                     K, ALPHA, V, LDV, T(K+1, 1), LDT)
                  CALL DTRMMOOP(SIDE, UPLO, 'No Transpose', TRANSV,
     $                     DIAGT, N-K, K, ALPHA, T(K+1, K+1), LDT,
     $                     V(K+1, 1), LDV, ONE, T(K+1, 1), LDT)
               END IF
            ELSE
               IF(VTRANS) THEN
                  CALL DTRMM('Left', 'Upper', TRANSV, DIAGV, N-K, K,
     $                     ALPHA, V(K+1, K+1), LDV, T(K+1, 1), LDT)
                  CALL DTRMMOOP(SIDE, UPLO, 'No Transpose', TRANSV,
     $                     DIAGT, N-K, K, ALPHA, T, LDT, V(1, K+1), LDV,
     $                     ONE, T(K+1, 1), LDT)
               ELSE
                  CALL DTRMM('Left', 'Lower', TRANSV, DIAGV, N-K, K,
     $                     ALPHA, V(K+1, K+1), LDV, T(K+1, 1), LDT)
                  CALL DTRMMOOP(SIDE, UPLO, 'No Transpose', TRANSV,
     $                     DIAGT, N-K, K, ALPHA, T, LDT, V(K+1, 1), LDV,
     $                     ONE, T(K+1, 1), LDT)
               END IF
            END IF
         END IF
         ! Since T_{11} and T_{22} are computed the same no matter what,
         ! we put the recursive calls here
         ! Compute T_{11} recursively
         CALL DTRTRM(SIDE, UPLO, TRANSV, DIAGT, DIAGV, K, ALPHA,
     $         T, LDT, V, LDV)
         ! Compute T_{22} recursively
         CALL DTRTRM(SIDE, UPLO, TRANSV, DIAGT, DIAGV, N-K, ALPHA,
     $         T(K+1, K+1), LDT, V(K+1, K+1), LDV)

      END SUBROUTINE
