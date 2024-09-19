c     Cost: n > k: 1/6 * (k^2-1)(2n+k)
c           n = k: 1/2 * (n^3-n)
      RECURSIVE SUBROUTINE MY_DLARFT_REC( DIRECT, STOREV, N, K, V, LDV,
     $                                    TAU, T, LDT)
         ! Arguemnts
         ! Scalars
         INTEGER           N, K, LDV, LDT
         CHARACTER         DIRECT, STOREV
         ! Matrix 
         DOUBLE PRECISION  V(LDV,*), T(LDT,*), TAU(N)

         ! Local variables
         INTEGER           I,J,L,INFO
         ! Parameters
         DOUBLE PRECISION ONE, NEG_ONE, ZERO
         PARAMETER(ONE=1.0D+0, ZERO = 0.0, NEG_ONE=-1.0D+0)

         ! Break V apart into 6 components
         ! V = |---------------|
         !     |V_{1,1} 0      |
         !     |V_{2,1} V_{2,2}|
         !     |V_{3,1} V_{3,2}|
         !     |---------------|
         ! V_{1,1}\in\R^{k,k} unit lower triangular
         ! V_{2,1}\in\R^{n-k,k} rectangular
         ! V_{3,1}\in\R^{m-n,k} rectangular
         ! 
         ! V_{2,2}\in\R^{n-k,n-k} unit upper triangular
         ! V_{3,2}\in\R^{m-n,n-k} rectangular

         ! We will construct the T matrix 
         ! T = |---------------| =  |--------|
         !     |T_{1,1} T_{1,2}|    |T_1  T_3|
         !     |0       T_{2,2}|    |0    T_2|
         !     |---------------|    |--------|

         ! T is the triangular factor attained from block reflectors. 
         ! To motivate the structure, consider the product
         !
         ! (I - V_1T_1V_1^\top)(I - V_2T_2V_2^\top)
         ! = I - V_1T_1V_1^\top - V_2T_2V_2^\top + V_1T_1V_1^\topV_2T_2V_2^\top
         !
         ! Define T_3 = -T_1V_1^\topV_2T_2
         !
         ! Then, we can define the matrix V as 
         ! V = |-------|
         !     |V_1 V_2|
         !     |-------|
         !
         ! So, our product is equivalent to the matrix product
         ! I - VTV^\top
         ! So, we compute T_1, then T_2, then use these values to get T_3
         !
         ! The general scheme used is inspired by the approach inside DGEQRT3
         ! which was (at the time of writing this code):
         ! Based on the algorithm of Elmroth and Gustavson,
         ! IBM J. Res. Develop. Vol 44 No. 4 July 2000.

         IF(K.EQ.0) THEN
            RETURN
         END IF
         ! Base case
         IF(K.EQ.1) THEN
            T(1,1) = TAU(1)
            RETURN
         END IF

         ! Beginning of executable statements
         L = K / 2

         ! Compute T_1
         CALL MY_DLARFT_REC(DIRECT, STOREV, N, L, V, LDV, TAU, T, LDT)

         ! Compute T_2
         CALL MY_DLARFT_REC(DIRECT, STOREV, N-L, K-L, V(L+1,L+1),
     $         LDV, TAU(L+1), T(L+1,L+1), LDT)

         ! Compute T_3 
         ! T_3 = V_{2,1}^\top
         DO J = 1, L
            DO I = L+1, K
               T(J,I) = V(I,J)
            END DO
         END DO
         ! T_3 = V_{2,1}^\top * V_{2,2}
         CALL DTRMM('Right', 'Lower', 'No transpose', 'Unit', 
     $         L, K - L, ONE, V(L+1, L+1), LDV, T(1, L+1), LDT)

         IF(N.GT.K) THEN
         ! T_3 = T_3 + V_{3,1}^\topV_{3,2}
            CALL DGEMM('Transpose', 'No transpose', L, K-L, N-K, ONE,
     $            V(K+1, 1), LDV, V(K+1,L+1), LDV, ONE, T(1, L+1), LDT)
         END IF

         ! At this point, we have that T_3 = V_1^\top *V_2
         ! All that is left is to pre and post multiply by -T_1 and T_2
         ! respectively.

         ! T_3 = -T_1*T_3
         CALL DTRMM('Left', 'Upper', 'No transpose', 'Non-unit',
     $         L, K - L, NEG_ONE, T, LDT, T(1, L+1), LDT)
         ! T_3 = T_3*T_2
         CALL DTRMM('Right', 'Upper', 'No transpose', 'Non-unit',
     $         L, K - L, ONE, T(L+1,L+1), LDT, T(1, L+1), LDT)
         
         ! Now, we have T in the correct form
      END SUBROUTINE
