c     Cost: 2/3 * (n^3 - n)
      RECURSIVE SUBROUTINE ZLUMM(N, ALPHA, Q, LDQ)
         ! Scalar Arguments
         INTEGER           N, LDQ
         COMPLEX*16  ALPHA
         ! Matrix Arguments
         COMPLEX*16  Q(LDQ,*)

         ! Local Variables
         ! Scalars
         INTEGER           K

         ! External Subroutines
         EXTERNAL          ZGEMM, ZTRMM

         ! Parameters
         COMPLEX*16  ONE
         PARAMETER(ONE=1.0D+0)

         ! Early exit if possible
         IF (N.EQ.0) THEN
            RETURN
         END IF


         ! Q is of the form
         !
         !     |---------------|
         ! Q = |Q_{1,1} Q_{1,2}|
         !     |Q_{2,1} Q_{2,2}|
         !     |---------------|
         ! k = n/2
         ! Q_{1,1}\in\R^{k\times k}
         ! Q_{1,2}\in\R^{k\times n-k}
         ! Q_{2,1}\in\R^{n-k\times k}
         ! Q_{2,2}\in\R^{n-k\times n-k}
         !
         ! Q is also of the form 
         !     |-|
         ! Q = |U|
         !     |L|
         !     |-|
         ! Due to the breakdown, Q_{1,1} and Q_{2,2} will also have the above
         ! structure
         ! Q_{1,2} is solely a part of U and Q_{2,1} is solely a part of L
         ! U is upper triangular and L is unit lower triangular

         ! U has the form
         ! U = |---------------|
         !     |U_{1,1} U_{1,2}|
         !     |0       U_{2,2}|
         !     |---------------|
         ! L has the form
         ! L = |---------------|
         !     |L_{1,1} 0      |
         !     |L_{2,1} L_{2,2}|
         !     |---------------|
         !
         ! So, we are computing
         ! Q = ALPHA*L*U = 
         !  |--------------------------------------------------------------------|
         !  |ALPHA*L_{1,1}*U_{1,1}  ALPHA*L_{1,1}*U_{1,2}                        |
         !  |ALPHA*L_{2,1}*U_{1,1}  ALPHA*L_{2,1}*U_{1,2} + ALPHA*L_{2,2}*U_{2,2}|
         !  |--------------------------------------------------------------------|
         ! IE
         !  Q_{1,1} = ALPHA*L_{1,1}*U_{1,1} (LUMM)
         !  Q_{1,2} = ALPHA*L_{1,1}*U_{1,2} (TRMM)
         !  Q_{2,1} = ALPHA*L_{2,1}*U_{1,1} (TRMM)
         !  Q_{2,2} = ALPHA*L_{2,1}*U_{1,2} + ALPHA*L_{2,2}*U_{2,2} (LUMM then GEMM)
         !  We compute these from bottom to top

         ! Base case of when N = 1
         IF (N.EQ.1) THEN
            ! We have a 1x1 matrix so we are multiplying a unit lower triangular
            ! matrix by an upper triangular matrix times a scalar. So we have
            ! that
            ! Q = ALPHA*L*U = ALPHA * U = ALPHA * Q
            Q(1,1) = ALPHA * Q(1,1)
            RETURN
         END IF
         K = N / 2

         ! Recursive Case
         ! Compute Q_{2,2} first
         ! Q_{2,2} = L_{2,2}*U_{2,2} (LUMM)
         CALL ZLUMM(N-K, ALPHA, Q(K + 1, K + 1), LDQ)

         ! Q_{2,2} = L_{2,1}*U_{1,2} + Q_{2,2} (GEMM)
         CALL ZGEMM('No transpose', 'No transpose', N-K, N-K, K, 
     $            ALPHA, Q(K+1,1), LDQ, Q(1,K+1), LDQ, ONE, 
     $            Q(K+1, K+1), LDQ)
         
         ! Compute Q_{2,1}
         ! Q_{2,1} = L_{2,1}*U_{1,1} (TRMM)
         CALL ZTRMM('Right', 'Upper', 'No-transpose', 'Non-unit',
     $            N-K, K, ALPHA, Q, LDQ, Q(K+1,1), LDQ)

         ! Compute Q_{2,1}
         CALL ZTRMM('Left', 'Lower', 'No-transpose', 'Unit',
     $            K, N-K, ALPHA, Q, LDQ, Q(1,K+1), LDQ)
         
         ! Compute Q_{1,1}
         CALL ZLUMM(K, ALPHA, Q, LDQ)

         ! Done!
      END SUBROUTINE
