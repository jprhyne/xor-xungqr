      RECURSIVE SUBROUTINE DTRMMOOP(SIDE, UPLO, TRANSA, TRANSB,
     $         DIAG, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
*
*        .. Scalar Arguments ..
         DOUBLE PRECISION  ALPHA, BETA
         INTEGER           M, N, LDA, LDB, LDC
         CHARACTER         SIDE, UPLO, TRANSB, TRANSC, DIAG
*        ..
*        .. Array Arguments ..
         DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*)
*        ..
*
*  =====================================================================
*
*        .. External Functions ..
         LOGICAL  LSAME
         EXTERNAL LSAME
*        ..
*        .. External Subroutines ..
*        ..
*        .. Intrinsic Functions ..
*        ..
*        .. Local Scalars ..
         INTEGER  L, K
         LOGICAL  LSIDE, UPPER, UNIT, TRANST, TRANSG
*        ..
*        .. Local Parameters ..
         DOUBLE PRECISION  ONE
         PARAMETER(ONE=1.0D+0)

*        ..
*
*        Beginning of Executable Statements
*
         LSIDE = LSAME(SIDE, 'L')
         UPPER = LSAME(UPLO, 'U')
         ! If we are transposing the triangular matrix (A)
         TRANST= LSAME(TRANSA, 'T').OR.LSAME(TRANSA, 'C')
         ! If we are transposing the general matrix (B)
         TRANSG= LSAME(TRANSB, 'T').OR.LSAME(TRANSB, 'C')
*
*        Terminating Case
*
         UNIT  = LSAME(DIAG, 'U')
*
*        Recursive Case
*
         L = MIN(M,N)
         K = L/2
         IF (LSIDE) THEN
*
*           We are multiplying A from the left IE we are computing
*           C = \alpha op(B)*op(A) + \beta C
*
            IF (UPPER) THEN
*
*              B is upper triangular
*
               IF (TRANST) THEN
*
*                 We are transposing B
*
                  IF (TRANSG) THEN
*
*                    We are transposing A.
*
*                    So we are computing
*                    C = \alpha B**T * A**T + \beta C. We break this down as follows
*
*                          |-------------|         |-------------------|
*                    C =   |C_{11} C_{12}| B**T =  |B_{11}**T 0        |
*                          |C_{21} C_{22}|         |B_{12}**T B_{22}**T|
*                          |-------------|         |-------------------|
*
*                          |-------------------|
*                    A**T =|A_{11}**T A_{21}**T|
*                          |A_{12}**T A_{22}**T|
*                          |-------------------|
*
*                    Which means that we get
*                    C_{11} = \alpha B_{11}**T * A_{11}**T + \beta C_{11}
*                    C_{12} = \alpha B_{11}**T * A_{21}**T + \beta C_{12}
*                    C_{21} = \alpha B_{12}**T * A_{11}**T + \alpha B_{22}**T * A_{12}**T + \beta C_{21}
*                    C_{22} = \alpha B_{12}**T * A_{21}**T + \alpha B_{22}**T * A_{22}**T + \beta C_{22}
*
*                    Computing C_{11} and C_{12} is just a recursive call to
*                    this routine but we can break down computing
*                    C_{21} and C_{22} as follows
*
*                    C_{21} = \alpha B_{12}**T * A_{11}**T + \beta C_{21} (GEMM call)
*                    C_{21} = \alpha B_{22}**T * A_{12}**T + C_{21} (This routine)
*
*                    C_{22} = \alpha B_{12}**T * A_{21}**T + \beta C_{22} (GEMM call)
*                    C_{22} = \alpha B_{22}**T * A_{22}**T + C_{22} (This routine)
*
                  ELSE
*
*                    We are not transposing A.
*
*                    So we are computing
*                    C = \alpha B**T * A + \beta C. We break this down as follows
*
*                          |-------------|         |-------------------|
*                    C =   |C_{11} C_{12}| B**T =  |B_{11}**T 0        |
*                          |C_{21} C_{22}|         |B_{12}**T B_{22}**T|
*                          |-------------|         |-------------------|
*
*                          |-------------|
*                    A =   |A_{11} A_{12}|
*                          |A_{21} A_{22}|
*                          |-------------|
*
*                    Which means that we get
*                    C_{11} = \alpha B_{11}**T * A_{11} + \beta C_{11}
*                    C_{12} = \alpha B_{11}**T * A_{12} + \beta C_{12}
*                    C_{21} = \alpha B_{12}**T * A_{11} + \alpha B_{22}**T * A_{21} + \beta C_{21}
*                    C_{22} = \alpha B_{12}**T * A_{12} + \alpha B_{22}**T * A_{22} + \beta C_{22}
*
*                    Computing C_{11} and C_{12} is just a recursive call to
*                    this routine but we can break down computing
*                    C_{21} and C_{22} as follows
*
*                    C_{21} = \alpha B_{12}**T * A_{11} + \beta C_{21} (GEMM call)
*                    C_{21} = \alpha B_{22}**T * A_{21} + C_{21} (This routine)
*
*                    C_{22} = \alpha B_{12}**T * A_{12} + \beta C_{22} (GEMM call)
*                    C_{22} = \alpha B_{22}**T * A_{22} + C_{22} (This routine)
*
                  ENDIF
               ELSE
*
*                 We are not transposing B
*
                  IF (TRANSG) THEN
*
*                    We are transposing A.
*
*                    So we are computing
*                    C = \alpha B * A**T + \beta C. We break this down as follows
*
*                          |-------------|      |-------------|
*                    C =   |C_{11} C_{12}| B =  |B_{11} B_{12}|
*                          |C_{21} C_{22}|      |0      B_{22}|
*                          |-------------|      |-------------|
*
*                          |-------------------|
*                    A**T =|A_{11}**T A_{21}**T|
*                          |A_{12}**T A_{22}**T|
*                          |-------------------|
*
*                    Which means that we get
*                    C_{11} = \alpha B_{11} * A_{11}**T + \alpha B_{12} * A_{12}**T + \beta C_{11}
*                    C_{12} = \alpha B_{11} * A_{21}**T + \alpha B_{12} * A_{22}**T + \beta C_{12}
*                    C_{21} = \alpha B_{22} * A_{12}**T + \beta C_{21}
*                    C_{22} = \alpha B_{22} * A_{22}**T + \beta C_{22}
*
*                    Computing C_{21} and C_{22} is just a recursive call to
*                    this routine but we can break down computing
*                    C_{11} and C_{12} as follows
*
*                    C_{11} = \alpha B_{12} * A_{12}**T + \beta C_{11} (GEMM call)
*                    C_{11} = \alpha B_{11} * A_{11}**T + C_{11} (This routine)
*
*                    C_{12} = \alpha B_{12} * A_{22}**T + \beta C_{12} (GEMM call)
*                    C_{12} = \alpha B_{11} * A_{21}**T + C_{12} (This routine)
*
                  ELSE
*
*                    We are not transposing A.
*
*                    So we are computing
*                    C = \alpha B * A + \beta C. We break this down as follows
*
*                          |-------------|      |-------------|
*                    C =   |C_{11} C_{12}| B =  |B_{11} B_{12}|
*                          |C_{21} C_{22}|      |0      B_{22}|
*                          |-------------|      |-------------|
*
*                          |-------------|
*                    A =   |A_{11} A_{12}|
*                          |A_{21} A_{22}|
*                          |-------------|
*
*                    Which means that we get
*                    C_{11} = \alpha B_{11} * A_{11} + \alpha B_{12} * A_{21} + \beta C_{11}
*                    C_{12} = \alpha B_{11} * A_{12} + \alpha B_{12} * A_{22} + \beta C_{12}
*                    C_{21} = \alpha B_{22} * A_{21} + \beta C_{21}
*                    C_{22} = \alpha B_{22} * A_{22} + \beta C_{22}
*
*                    Computing C_{21} and C_{22} is just a recursive call to
*                    this routine but we can break down computing
*                    C_{11} and C_{12} as follows
*
*                    C_{11} = \alpha B_{12} * A_{21} + \beta C_{11} (GEMM call)
*                    C_{11} = \alpha B_{11} * A_{11} + C_{11} (This routine)
*
*                    C_{12} = \alpha B_{12} * A_{22} + \beta C_{12} (GEMM call)
*                    C_{12} = \alpha B_{11} * A_{12} + C_{12} (This routine)
*
                  ENDIF
               END IF
            ELSE
*
*              B is lower triangular
*
               IF (TRANST) THEN
*
*                 We are transposing B
*
                  IF (TRANSG) THEN
*
*                    We are transposing A.
*
*                    So we are computing
*                    C = \alpha B**T * A**T + \beta C. We break this down as follows
*
*                          |-------------|         |-------------------|
*                    C =   |C_{11} C_{12}| B**T =  |B_{11}**T B_{21}**T|
*                          |C_{21} C_{22}|         |0         B_{22}**T|
*                          |-------------|         |-------------------|
*
*                          |-------------------|
*                    A**T =|A_{11}**T A_{21}**T|
*                          |A_{12}**T A_{22}**T|
*                          |-------------------|
*
*                    Which means that we get
*                    C_{11} = \alpha B_{11}**T * A_{11}**T + \alpha B_{21}**T * A_{12}**T + \beta C_{11}
*                    C_{12} = \alpha B_{11}**T * A_{21}**T + \alpha B_{21}**T * A_{22}**T + \beta C_{12}
*                    C_{21} = \alpha B_{22}**T * A_{12}**T + \beta C_{21}
*                    C_{22} = \alpha B_{22}**T * A_{22}**T + \beta C_{22}
*
*                    Computing C_{21} and C_{22} is just a recursive call to
*                    this routine but we can break down computing
*                    C_{11} and C_{12} as follows
*
*                    C_{11} = \alpha B_{21}**T * A_{12}**T + \beta C_{11} (GEMM call)
*                    C_{11} = \alpha B_{11}**T * A_{11}**T + C_{11} (This routine)
*
*                    C_{12} = \alpha B_{21}**T * A_{22}**T + \beta C_{12} (GEMM call)
*                    C_{12} = \alpha B_{11}**T * A_{21}**T + C_{12} (This routine)
*
                  ELSE
*
*                    We are not transposing A.
*
*                    So we are computing
*                    C = \alpha B**T * A + \beta C. We break this down as follows
*
*                          |-------------|         |-------------------|
*                    C =   |C_{11} C_{12}| B**T =  |B_{11}**T B_{21}**T|
*                          |C_{21} C_{22}|         |0         B_{22}**T|
*                          |-------------|         |-------------------|
*
*                          |-------------|
*                    A =   |A_{11} A_{12}|
*                          |A_{21} A_{22}|
*                          |-------------|
*
*                    Which means that we get
*                    C_{11} = \alpha B_{11}**T * A_{11} + \alpha B_{21}**T * A_{21} + \beta C_{11}
*                    C_{12} = \alpha B_{11}**T * A_{12} + \alpha B_{21}**T * A_{22} + \beta C_{12}
*                    C_{21} = \alpha B_{22}**T * A_{21} + \beta C_{21}
*                    C_{22} = \alpha B_{22}**T * A_{22} + \beta C_{22}
*
*                    Computing C_{21} and C_{22} is just a recursive call to
*                    this routine but we can break down computing
*                    C_{11} and C_{12} as follows
*
*                    C_{11} = \alpha B_{21}**T * A_{21} + \beta C_{11} (GEMM call)
*                    C_{11} = \alpha B_{11}**T * A_{11} + C_{11} (This routine)
*
*                    C_{12} = \alpha B_{21}**T * A_{22} + \beta C_{12} (GEMM call)
*                    C_{12} = \alpha B_{11}**T * A_{12} + C_{12} (This routine)
*
                  ENDIF
               ELSE
*
*                 We are not transposing B
*
                  IF (TRANSG) THEN
*
*                    We are transposing A.
*
*                    So we are computing
*                    C = \alpha B * A**T + \beta C. We break this down as follows
*
*                          |-------------|      |-------------|
*                    C =   |C_{11} C_{12}| B =  |B_{11} 0     |
*                          |C_{21} C_{22}|      |B_{21} B_{22}|
*                          |-------------|      |-------------|
*
*                          |-------------------|
*                    A**T =|A_{11}**T A_{21}**T|
*                          |A_{12}**T A_{22}**T|
*                          |-------------------|
*
*                    Which means that we get
*                    C_{11} = \alpha B_{11} * A_{11}**T + \beta C_{11}
*                    C_{12} = \alpha B_{11} * A_{21}**T + \beta C_{12}
*                    C_{21} = \alpha B_{21} * A_{11}**T + \alpha B_{22} * A_{12}**T + \beta * C_{21}
*                    C_{22} = \alpha B_{21} * A_{21}**T + \alpha B_{22} * A_{22}**T + \beta * C_{22}
*
*                    Computing C_{11} and C_{12} is just a recursive call to
*                    this routine but we can break down computing
*                    C_{21} and C_{22} as follows
*
*                    C_{21} = \alpha B_{21} * A_{11}**T + \beta C_{21} (GEMM call)
*                    C_{21} = \alpha B_{22} * A_{12}**T + C_{21} (This routine)
*
*                    C_{22} = \alpha B_{21} * A_{21}**T + \beta C_{22} (GEMM call)
*                    C_{22} = \alpha B_{22} * A_{22}**T + C_{22} (This routine)
*
                  ELSE
*
*                    We are not transposing A.
*
*                    So we are computing
*                    C = \alpha B * A + \beta C. We break this down as follows
*
*                          |-------------|      |-------------|
*                    C =   |C_{11} C_{12}| B =  |B_{11} 0     |
*                          |C_{21} C_{22}|      |B_{21} B_{22}|
*                          |-------------|      |-------------|
*
*                          |-------------|
*                    A =   |A_{11} A_{12}|
*                          |A_{21} A_{22}|
*                          |-------------|
*
*                    Which means that we get
*                    C_{11} = \alpha B_{11} * A_{11} + \beta C_{11}
*                    C_{12} = \alpha B_{11} * A_{12} + \beta C_{12}
*                    C_{21} = \alpha B_{21} * A_{11} + \alpha B_{22} * A_{21} + \beta * C_{21}
*                    C_{22} = \alpha B_{21} * A_{12} + \alpha B_{22} * A_{22} + \beta * C_{22}
*
*                    Computing C_{11} and C_{12} is just a recursive call to
*                    this routine but we can break down computing
*                    C_{21} and C_{22} as follows
*
*                    C_{21} = \alpha B_{21} * A_{11} + \beta C_{21} (GEMM call)
*                    C_{21} = \alpha B_{22} * A_{21} + C_{21} (This routine)
*
*                    C_{22} = \alpha B_{21} * A_{12} + \beta C_{22} (GEMM call)
*                    C_{22} = \alpha B_{22} * A_{22} + C_{22} (This routine)
*
                  ENDIF
               END IF
            END IF
         ELSE
*
*           We are multiplying A from the right IE we are computing
*           C = \alpha op(A)*op(B) + \beta C
*
            IF (UPPER) THEN
*
*              B is upper triangular
*
               IF (TRANST) THEN
*
*                 We are transposing B
*
                  IF (TRANSG) THEN
*
*                    We are transposing A.
*
*                    So we are computing
*                    C = \alpha  A**T * B**T + \beta C. We break this down as follows
*
*                          |-------------|         |-------------------|
*                    C =   |C_{11} C_{12}| B**T =  |B_{11}**T 0        |
*                          |C_{21} C_{22}|         |B_{12}**T B_{22}**T|
*                          |-------------|         |-------------------|
*
*                          |-------------------|
*                    A**T =|A_{11}**T A_{21}**T|
*                          |A_{12}**T A_{22}**T|
*                          |-------------------|
*
*                    Which means that we get
*                    C_{11} = \alpha A_{11}**T * B_{11}**T + \alpha A_{21}**T * B_{12}**T + \beta C_{11}
*                    C_{12} = \alpha A_{21}**T * B_{22}**T + \beta C_{12}
*                    C_{21} = \alpha A_{12}**T * B_{11}**T + \alpha A_{22}**T * B_{12}**T + \beta C_{21}
*                    C_{22} = \alpha A_{22}**T * B_{22}**T + \beta C_{22}
*
*                    Computing C_{12} and C_{22} is just a recursive call to
*                    this routine but we can break down computing
*                    C_{11} and C_{21} as follows
*
*                    C_{11} = \alpha A_{21}**T * B_{12}**T + \beta C_{11} (GEMM call)
*                    C_{11} = \alpha A_{11}**T * B_{11}**T + C_{11} (This routine)
*
*                    C_{21} = \alpha A_{22}**T * B_{12}**T + \beta C_{21} (GEMM call)
*                    C_{21} = \alpha A_{12}**T * B_{11}**T + C_{21} (This routine)
*
                  ELSE
*
*                    We are not transposing A.
*
*                    So we are computing
*                    C = \alpha A * B**T + \beta C. We break this down as follows
*
*                          |-------------|         |-------------------|
*                    C =   |C_{11} C_{12}| B**T =  |B_{11}**T 0        |
*                          |C_{21} C_{22}|         |B_{12}**T B_{22}**T|
*                          |-------------|         |-------------------|
*
*                          |-------------|
*                    A =   |A_{11} A_{12}|
*                          |A_{21} A_{22}|
*                          |-------------|
*
*                    Which means that we get
*                    C_{11} = \alpha A_{11} * B_{11}**T + \alpha A_{12} * B_{12}**T + \beta C_{11}
*                    C_{12} = \alpha A_{12} * B_{22}**T + \beta C_{12}
*                    C_{21} = \alpha A_{21} * B_{11}**T + \alpha A_{22} * B_{12}**T + \beta C_{21}
*                    C_{22} = \alpha A_{22} * B_{22}**T + \beta C_{22}
*
*                    Computing C_{12} and C_{22} is just a recursive call to
*                    this routine but we can break down computing
*                    C_{11} and C_{21} as follows
*
*                    C_{11} = \alpha A_{12} * B_{12}**T + \beta C_{11} (GEMM call)
*                    C_{11} = \alpha A_{11} * B_{11}**T + C_{11} (This routine)
*
*                    C_{21} = \alpha A_{22} * B_{12}**T + \beta C_{21} (GEMM call)
*                    C_{21} = \alpha A_{21} * B_{11}**T + C_{21} (This routine)
*
                  ENDIF
               ELSE
*
*                 We are not transposing B
*
                  IF (TRANSG) THEN
*
*                    We are transposing A.
*
*                    So we are computing
*                    C = \alpha A**T * B + \beta C. We break this down as follows
*
*                          |-------------|      |-------------|
*                    C =   |C_{11} C_{12}| B =  |B_{11} B_{12}|
*                          |C_{21} C_{22}|      |0      B_{22}|
*                          |-------------|      |-------------|
*
*                          |-------------------|
*                    A**T =|A_{11}**T A_{21}**T|
*                          |A_{12}**T A_{22}**T|
*                          |-------------------|
*
*                    Which means that we get
*                    C_{11} = \alpha A_{11}**T * B_{11} + \beta C_{11}
*                    C_{12} = \alpha A_{11}**T * B_{12} + \alpha A_{21}**T * B_{22} + \beta C_{12}
*                    C_{21} = \alpha A_{12}**T * B_{11} + \beta C_{21}
*                    C_{22} = \alpha A_{12}**T * B_{12} + \alpha A_{22}**T * B_{22} + \beta C_{22}
*
*                    Computing C_{11} and C_{21} is just a recursive call to
*                    this routine but we can break down computing
*                    C_{12} and C_{22} as follows
*
*                    C_{12} = \alpha A_{11}**T * B_{12} + \beta C_{12} (GEMM call)
*                    C_{12} = \alpha A_{21}**T * B_{22} + C_{12} (This routine)
*
*                    C_{22} = \alpha A_{12}**T * B_{12} + \beta C_{22} (GEMM call)
*                    C_{22} = \alpha A_{22}**T * B_{22} + C_{22} (This routine)
*
                  ELSE
*
*                    We are not transposing A.
*
*                    So we are computing
*                    C = \alpha A * B + \beta C. We break this down as follows
*
*                          |-------------|      |-------------|
*                    C =   |C_{11} C_{12}| B =  |B_{11} B_{12}|
*                          |C_{21} C_{22}|      |0      B_{22}|
*                          |-------------|      |-------------|
*
*                          |-------------|
*                    A =   |A_{11} A_{12}|
*                          |A_{21} A_{22}|
*                          |-------------|
*
*                    Which means that we get
*                    C_{11} = \alpha A_{11} * B_{11} + \beta C_{11}
*                    C_{12} = \alpha A_{11} * B_{12} + \alpha A_{12} * B_{22} + \beta C_{12}
*                    C_{21} = \alpha A_{21} * B_{11} + \beta C_{21}
*                    C_{22} = \alpha A_{21} * B_{12} + \alpha A_{22} * B_{22} + \beta C_{22}
*
*                    Computing C_{11} and C_{21} is just a recursive call to
*                    this routine but we can break down computing
*                    C_{12} and C_{22} as follows
*
*                    C_{12} = \alpha A_{11} * B_{12} + \beta C_{12} (GEMM call)
*                    C_{12} = \alpha A_{12} * B_{22} + C_{12} (This routine)
*
*                    C_{22} = \alpha A_{21} * B_{12} + \beta C_{22} (GEMM call)
*                    C_{22} = \alpha A_{22} * B_{22} + C_{22} (This routine)
*
                  ENDIF
               END IF
            ELSE
*
*              B is lower triangular
*
               IF (TRANST) THEN
*
*                 We are transposing B
*
                  IF (TRANSG) THEN
*
*                    We are transposing A.
*
*                    So we are computing
*                    C = \alpha A**T * B**T + \beta C. We break this down as follows
*
*                          |-------------|         |-------------------|
*                    C =   |C_{11} C_{12}| B**T =  |B_{11}**T B_{21}**T|
*                          |C_{21} C_{22}|         |0         B_{22}**T|
*                          |-------------|         |-------------------|
*
*                          |-------------------|
*                    A**T =|A_{11}**T A_{21}**T|
*                          |A_{12}**T A_{22}**T|
*                          |-------------------|
*
*                    Which means that we get
*                    C_{11} = \alpha A_{11}**T * B_{11} + \beta C_{11}
*                    C_{12} = \alpha A_{11}**T * B_{21}**T + \alpha A_{21}**T * B_{22}**T + \beta C_{12}
*                    C_{21} = \alpha A_{12}**T * B_{11}**T + \beta C_{21}
*                    C_{22} = \alpha A_{12}**T * B_{21}**T + \alpha A_{22}**T * B_{22}**T + \beta C_{22}
*
*                    Computing C_{11} and C_{21} is just a recursive call to
*                    this routine but we can break down computing
*                    C_{12} and C_{22} as follows
*
*                    C_{12} = \alpha A_{11}**T * B_{21}**T + \beta C_{12} (GEMM call)
*                    C_{12} = \alpha A_{21}**T * B_{22}**T + C_{12} (This routine)
*
*                    C_{22} = \alpha A_{12}**T * B_{21}**T + \beta C_{22} (GEMM call)
*                    C_{22} = \alpha A_{22}**T * B_{22}**T + C_{22} (This routine)
*
                  ELSE
*
*                    We are not transposing A.
*
*                    So we are computing
*                    C = \alpha A * B**T + \beta C. We break this down as follows
*
*                          |-------------|         |-------------------|
*                    C =   |C_{11} C_{12}| B**T =  |B_{11}**T B_{21}**T|
*                          |C_{21} C_{22}|         |0         B_{22}**T|
*                          |-------------|         |-------------------|
*
*                          |-------------|
*                    A =   |A_{11} A_{12}|
*                          |A_{21} A_{22}|
*                          |-------------|
*
*                    Which means that we get
*                    C_{11} = \alpha A_{11} * B_{11} + \beta C_{11}
*                    C_{12} = \alpha A_{11} * B_{21}**T + \alpha A_{12} * B_{22}**T + \beta C_{12}
*                    C_{21} = \alpha A_{21} * B_{11}**T + \beta C_{21}
*                    C_{22} = \alpha A_{21} * B_{21}**T + \alpha A_{22} * B_{22}**T + \beta C_{22}
*
*                    Computing C_{11} and C_{21} is just a recursive call to
*                    this routine but we can break down computing
*                    C_{12} and C_{22} as follows
*
*                    C_{12} = \alpha A_{11} * B_{21}**T + \beta C_{12} (GEMM call)
*                    C_{12} = \alpha A_{12} * B_{22}**T + C_{12} (This routine)
*
*                    C_{22} = \alpha A_{21} * B_{21}**T + \beta C_{22} (GEMM call)
*                    C_{22} = \alpha A_{22} * B_{22}**T + C_{22} (This routine)
*
                  ENDIF
               ELSE
*
*                 We are not transposing B
*
                  IF (TRANSG) THEN
*
*                    We are transposing A.
*
*                    So we are computing
*                    C = \alpha A**T * B + \beta C. We break this down as follows
*
*                          |-------------|      |-------------|
*                    C =   |C_{11} C_{12}| B =  |B_{11} 0     |
*                          |C_{21} C_{22}|      |B_{21} B_{22}|
*                          |-------------|      |-------------|
*
*                          |-------------------|
*                    A**T =|A_{11}**T A_{21}**T|
*                          |A_{12}**T A_{22}**T|
*                          |-------------------|
*
*                    Which means that we get
*                    C_{11} = \alpha A_{11}**T * B_{11} + \alpha A_{21}**T * B_{21} + \beta C_{11}
*                    C_{12} = \alpha A_{21}**T * B_{22} + \beta C_{12}
*                    C_{21} = \alpha A_{12}**T * B_{11} + \alpha A_{22}**T * B_{21} + \beta C_{21}
*                    C_{22} = \alpha A_{22}**T * B_{22} + \beta C_{22}
*
*                    Computing C_{12} and C_{22} is just a recursive call to
*                    this routine but we can break down computing
*                    C_{11} and C_{21} as follows
*
*                    C_{11} = \alpha A_{21}**T * B_{21} + \beta C_{11} (GEMM call)
*                    C_{11} = \alpha A_{11}**T * B_{11} + C_{11}(This routine)
*
*                    C_{21} = \alpha A_{22}**T * B_{21} + \beta C_{21} (GEMM call)
*                    C_{21} = \alpha A_{12}**T * B_{11} + C_{21} (This routine)
*
                  ELSE
*
*                    We are not transposing A.
*
*                    So we are computing
*                    C = \alpha A * B + \beta C. We break this down as follows
*
*                          |-------------|      |-------------|
*                    C =   |C_{11} C_{12}| B =  |B_{11} 0     |
*                          |C_{21} C_{22}|      |B_{21} B_{22}|
*                          |-------------|      |-------------|
*
*                          |-------------|
*                    A =   |A_{11} A_{12}|
*                          |A_{21} A_{22}|
*                          |-------------|
*
*                    Which means that we get
*                    C_{11} = \alpha A_{11} * B_{11} + \alpha A_{12} * B_{21} + \beta C_{11}
*                    C_{12} = \alpha A_{12} * B_{22} + \beta C_{12}
*                    C_{21} = \alpha A_{21} * B_{11} + \alpha A_{22} * B_{21} + \beta C_{21}
*                    C_{22} = \alpha A_{22} * B_{22} + \beta C_{22}
*
*                    Computing C_{12} and C_{22} is just a recursive call to
*                    this routine but we can break down computing
*                    C_{11} and C_{21} as follows
*
*                    C_{11} = \alpha A_{12} * B_{21} + \beta C_{11} (GEMM call)
*                    C_{11} = \alpha A_{11} * B_{11} + C_{11}(This routine)
*
*                    C_{21} = \alpha A_{22} * B_{21} + \beta C_{21} (GEMM call)
*                    C_{21} = \alpha A_{21} * B_{11} + C_{21} (This routine)
*
                  ENDIF
               END IF
            END IF
         END IF
      END SUBROUTINE
