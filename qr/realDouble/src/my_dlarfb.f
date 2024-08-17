      SUBROUTINE MY_DLARFB(M, N, K, A, LDA, C, LDC)
         ! Scalar arguments
         INTEGER           M, N, K, LDA, LDC
         ! Array arguments
         DOUBLE PRECISION  A(LDA,*), C(LDC,*)
         ! Parameters
         DOUBLE PRECISION ONE, ZERO
         PARAMETER(ONE=1.0E+0, ZERO = 0.0E+0)
         ! On input, A will be of the form   A =   ( T )
         !                                         ( V )
         ! Where T is upper triangular, and V is broken further down into
         !                                   V =   ( V1 )
         !                                         ( V2 )
         ! Where V1 is unit lower triangular (diagonal assumed) and V2 is rectangular
         ! V1 and T start at A(1,1)
*
*           Apply H to A(i:m,i+ib:n) from the left
*
*
*              Form  H * C  or  H**T * C  where  C = ( C1=0 )
*                                                    ( C2=* )
*
*              C := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
*
*              C1 := V2**T * C2
*
         CALL DGEMM( 'Transpose', 'No transpose', K, N, M - K, ONE, 
     $                  A(1 + K, 1), LDA, C(1 + K, 1), LDC, ZERO, 
     $                  C, LDC)
*               CALL DGEMM( 'Transpose', 'No transpose', IB, 
*     $                        N-I-IB+1, M-I+1-IB,
*     $                        ONE, A( I+IB, I ), LDA, 
*     $                        A( I + IB, I + IB ), LDA,
*     $                        ZERO, A(I,I+IB), LDA )
*
**              W  := W * T**T  or  W * T
*              C1 := T * C1
*
         CALL DTRMM( 'Left', 'Upper', 'No transpose', 'Non-unit', K, N,
     $                  ONE, A, LDA, C, LDC)
*               CALL DTRMM( 'Left', 'Upper', 'No transpose', 'Non-unit',
*     $                     IB, N-I-IB+1, ONE, A(I,I), LDA, 
*     $                     A(I,I+IB), LDA )
*
*              C := C - V * W**T
*
*
*                 C2 := C2 - V2 * C1
*
         CALL DGEMM( 'No transpose', 'No transpose', M - K, N, K, -ONE,
     $                  A(1 + K, 1), LDA, C, LDC, ONE, C(1 + K, 1), LDC)
*                  CALL DGEMM( 'No transpose', 'No transpose', M-I-IB+1,
*     $                        N-I-IB+1, IB,
*     $                        -ONE, A( I+IB, I ), LDA, A(I,I+IB), LDA,
*     $                        ONE, A( I+IB, I+IB ), LDA )
*
*              C1 := -V1 * C1
*
         CALL DTRMM( 'Left', 'Lower', 'No transpose', 'Unit', K, N,
     $                  -ONE, A, LDA, C, LDC)
*               CALL DTRMM( 'Left', 'Lower', 'Non transpose', 'Unit', IB,
*     $                     N-I-IB+1, -ONE, A(I,I), LDA, A(I,I+IB), LDA )


      END SUBROUTINE
