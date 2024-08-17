* We compute C = T*A**T + ALPHA*C or C = A**T*T + ALPHA * C 
* If SIDE = 'L'/'l' or 'R'/'r' respectively
*     For SIDE = 'L', T is upper triangular and for SIDE = 'R', T is lower
*     triangular. This is due to where we need this functionality in our
*     codebase
* Currently do not support any other functionality, but can if desired
*     C is m by n
*     T is m by m
*     A is n by m -> A**T is m by n
      RECURSIVE SUBROUTINE ZTRMMOOP(SIDE, DIAG, M, N, A, LDA, T, LDT,
     $                              ALPHA, C, LDC, INFO)
*
*        .. Scalar Arguments ..
*
         INTEGER           M, N, LDA, LDT, LDC, INFO
         CHARACTER         SIDE, DIAG
         COMPLEX*16        ALPHA

*
*        .. Array arguments ..
*
         COMPLEX*16        A(LDA,*), T(LDT,*), C(LDC,*)
*
*        .. Local variables ..
*
         INTEGER           I,J,K
         COMPLEX*16        SUM, SCAL
*
*        .. Local parameters ..
*
         COMPLEX*16        ONE, ZERO
         PARAMETER(ONE=1.0d+0,ZERO=0.0d+0)
*
*        .. External Functions ..
*
         COMPLEX*16        ZDOTC
         EXTERNAL          ZDOTC
*
*        .. External Subroutines ..
*
         EXTERNAL          ZAXPY
*        We break down each matrix into the following form
*
*        |-------|   |-------|   |-------|   |-------|**T
*        |C11 C12| = |C11 C12| + |T11 T12| * |A11 A12|
*        |C21 C22|   |C21 C22|   |T21 T22|   |A21 A22|
*        |-------|   |-------|   |-------|   |-------|
*
*
*
*        C_{1,1}\in\R^{k by k}
*        C_{1,2}\in\R^{k by n - k}
*        C_{2,1}\in\R^{m - k by k}
*        C_{2,2}\in\R^{m - k by n - k}
*        
*        We choose K = MIN(M,N)/2
*
*        C_{2,2} = C_{2,2} + T_{2,2}*A_{2,2}**T
*        C_{2,1} = C_{2,1} + T_{2,2}*A_{1,2}**T
*
*        C_{1,1} = C_{1,1} + T_{1,1}*A_{1,1}**T + T_{1,2}*A_{1,2}**T
*        C_{1,2} = C_{1,2} + T_{1,1}*A_{2,1}**T + T_{1,2}*A_{2,2}**T
*
*        C_{1,1} can be broken down into two different operations
*        C_{1,1} = C_{1,1} + T_{1,2}*A_{1,2}**T [GEMM]
*        C_{1,1} = C_{1,1} + T_{1,1}*A_{1,1}**T [DTRMMOOP]
*
*        C_{1,2} can be broken down into two different operations
*        C_{1,2} = C_{1,2} + T_{1,2}*A_{2,2}**T [GEMM]
*        C_{1,2} = C_{1,2} + T_{1,1}*A_{2,1}**T [DTRMMOOP]

*--------------------------------------------------------------------------
*        Begin of executable statements
*--------------------------------------------------------------------------
*        Base case for the left half of the matrix
         IF (M.EQ.0.OR.N.EQ.0) THEN
            RETURN
         END IF
         ! Setting the return value now. Update to -1 iff side is invalid
         INFO = 0
*        Determine if we have T on the left or right
         IF(SIDE.EQ.'L'.OR.SIDE.EQ.'l') GOTO 10
         IF(SIDE.EQ.'R'.OR.SIDE.EQ.'r') GOTO 20
         INFO = -1
         RETURN
*        Base cases
   10    IF (M.EQ.1) THEN
*           In this case, T is 1x1 upper triangular matrix.
*           Therefore, we need to compute C = C + A*T(1,1)
*
*           This special case is done because for some reason, when we go to the
*           10 do loop, we have j=1,1 and this somehow goes to j=2. Not sure
*           why, but this is a workaround. (maybe some gdb issue. will toy
*           around with removing later as a last cleanup step).
*
            ! CALL DAXPY(N, T(1,1), A(1,1), 1, C(1,1), LDC)
            CALL ZSCAL(N, ALPHA, C(1,1), LDC)
            CALL ZAXPY(N, T(1,1),A(1,1), 1, C(1,1), LDC)
            RETURN
         ELSE IF (N.EQ.1) THEN
            ! Write quick implementation of dtrmvoop (should be similar to this
            ! case)
*           In this case, we have a C as a column vector and we need to compute
*           a modified matrix vector product (ie a modified form of dtrmv) But
*           also out of place.
*           We accomplish this by computing each element of C through ddot.
            
            DO I = 1,M
               C(I,1) = ALPHA * C(I,1) + ZDOTC(M-I+1, T(I,I), LDT, 
     $                     A(1,I),LDA)
            END DO
            RETURN
         END IF
*        Recursive case
         K = MIN(M,N) / 2
*        Compute C21
         CALL ZTRMMOOP(SIDE, DIAG, M - K, K, A(1, K+1), LDA,
     $            T(K + 1, K + 1), LDT, ALPHA, C(K + 1, 1), LDC, INFO)
*        Compute C22
         CALL ZTRMMOOP(SIDE, DIAG, M - K, N - K, A(K+1,K+1), LDA,
     $           T(K+1,K+1), LDT, ALPHA, C(K + 1, K + 1), LDC, INFO)
*        Compute C11 part 1
         CALL ZTRMMOOP(SIDE, DIAG, K, K, A, LDA, T, LDT, ALPHA, C, LDC, 
     $                  INFO)
*        Compute C11 part 2
         CALL ZGEMM('No transpose', 'Transpose', K, K, M - K, ONE, 
     $           T(1, K + 1), LDT, A(1, K + 1), LDA, ONE, C, LDC)
*        Compute C12 part 1
         CALL ZTRMMOOP(SIDE, DIAG, K, N-K, A(K + 1, 1), LDA, T, LDT, 
     $           ALPHA, C(1, K + 1), LDC, INFO)
*        Compute C12 part 2
         CALL ZGEMM('No transpose', 'Transpose', K, N - K, M - K,
     $           ONE, T(1, K + 1), LDT, A(K + 1, K + 1), LDA,
     $           ONE, C(1, K + 1), LDC)
         INFO = 0
         RETURN
   20    IF (N.EQ.1) THEN 
            CALL ZSCAL(M, ALPHA, C, 1)
            IF (DIAG.EQ.'U'.OR.DIAG.EQ.'u') THEN
               SCAL = ONE
            ELSE
               SCAL = T(1,1)
            END IF
            CALL ZAXPY(M, SCAL, A, LDA, C, 1)
            RETURN
         END IF
         IF (M.EQ.1) THEN
            CALL ZSCAL(N, ALPHA, C, LDC)
            ! N >= 2
            IF (DIAG.EQ.'U'.OR.DIAG.EQ.'u') THEN
               DO J = N, 1, -1
                  C(1,J) = C(1,J) + A(J,1) + ZDOTC(N - J, A(J + 1, 1),
     $               1, T(J + 1,J), 1)
               END DO
            ELSE
               DO J = N, 1, -1
                  C(1,J) = C(1,J)  + ZDOTC(N - J + 1, A(J, 1), 1, 
     $               T(J,J), 1)
               END DO
            END IF
            RETURN
         END IF
         K = MIN(M,N) / 2
         ! Compute C12
         CALL ZTRMMOOP(SIDE, DIAG, K, N-K, A(K+1, 1), LDA, T(K+1,K+1),
     $            LDT, ALPHA, C(1, K + 1), LDC, INFO)
         ! Compute C22
         CALL ZTRMMOOP(SIDE, DIAG, M-K, N-K, A(K+1,K+1), LDA,
     $            T(K+1,K+1), LDT, ALPHA, C(K+1,K+1), LDC, INFO)
         ! Compute C11 part 1
         CALL ZGEMM('Transpose', 'Non-transpose', K, K, N-K, ONE, 
     $            A(K+1,1), LDA, T(K+1,1), LDT, ALPHA, C, LDC)
         ! Compute C11 part 2
         CALL ZTRMMOOP(SIDE, DIAG, K, K, A, LDA, T, LDT, ONE, C, 
     $            LDC, INFO)
         ! Compute C21 part 1
         CALL ZGEMM('Transpose', 'Non-transpose', M-K, K, N-K, ONE,
     $            A(K+1,K+1), LDA, T(K+1, 1), LDT, ALPHA, C(K+1,1),LDC)
         ! Compute C21 part 2
         CALL ZTRMMOOP(SIDE, DIAG, M-K, K, A(1, K+1), LDA, T, LDT,
     $            ONE, C(K+1,1), LDC, INFO)
         RETURN
      END SUBROUTINE
