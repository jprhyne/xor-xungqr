*     This is a recursive subroutine that will compute T*V**T as part of
*     the dorg2r algorithm.
*
*     Parameters
*     N(in):     Number of columns in the matrix Q
*     Q(in/out): on input Matrix that will hold V and T as described below.
*        V is input only
*        T is input/output
*
*         |-|
*     Q = |T|
*         |-|
*         |V|
*         |-|
*
*     V = |---|
*         |V_1|
*         |---|
*         |V_2|
*         |---|
*        On output, we replace T with T*V_1**T
*
*
*
*     Cost: (2n**3 + 3n**2 -5n) / 6

      RECURSIVE SUBROUTINE ZTVT(N, Q, LDQ)
*
*        Input paramaters
*
         INTEGER           N, LDQ

         COMPLEX*16        Q(LDQ,*)
*
*        Local variables
*
         INTEGER           K, INFO
*
*        External subroutines
*
         EXTERNAL ZTRMM, ZTRMMOOP

*
*        Local parameters
*
         COMPLEX*16           ONE
         PARAMETER(ONE=1.0d+0)
*
*        Beginning of executable statements
*
*
*        Base case
*
         IF (N.LT.2) THEN
            RETURN
         ELSE
*
*        Recursive case
*
            K = N / 2
*           Computes T_{1,2} = T_{1,2}V_{2,2}^\top
            CALL ZTRMM('Right', 'Lower', 'Transpose', 'Unit', K, N - K, 
     $         ONE, Q(K + 1, K + 1), LDQ, Q(1, K + 1), LDQ)
*           Compute T_{1,2} = T_{1,2} + T_{1,1}V_{2,1}^\top
            CALL ZTRMMOOP('Left', 'N', K, N - K, Q(K + 1, 1), LDQ, Q,
     $         LDQ, ONE, Q(1, K + 1), LDQ, INFO)

*           Compute T_{1,1} = T_{1,1}V_{1,1}^\top
            CALL ZTVT(K, Q, LDQ)

*           Compute T_{2,2} = T_{2,2}V_{2,2}^\top
            CALL ZTVT(N-K, Q(K + 1, K + 1), LDQ)
         END IF

            
      END SUBROUTINE 
