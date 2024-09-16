c     Cost: m > n: 1/6 * (n^2-1)(2m+n)
c           m = n: 1/2 * (n^3-n)
      RECURSIVE SUBROUTINE MY_DLARFT_REC(DIRECT, STOREV, M, N, V, LDV, 
     $                                   TAU, T, LDT)
         ! Arguments
         ! Scalars
         CHARACTER         DIRECT, STOREV
         INTEGER           M, N, LDV, LDT
         ! Matrix 
         DOUBLE PRECISION  V(LDV,*), T(LDT,*), TAU(N)
         ! External subroutines
         EXTERNAL          DGEMM, DTRMM
         ! External functions
         LOGICAL           LSAME
         EXTERNAL          LSAME

         ! Local variables
         INTEGER           I,J,K,INFO,V1I,V1J,T3I,T3J,AI,AJ,BI,BJ
         INTEGER           TLEFTI,TLEFTJ,TRIGHTI,TRIGHTJ
         LOGICAL           COLV, DIRF, COPYTR
         ! TRMMS is the value of side that is used in the first trmm call
         ! TRMMT is the value of trans that is used in the frist trmm call
         ! TRMM1U is the value of uplo that is used in the first trmm call
         ! TRMM2U is the value of uplo that is used in the second trmm call
         ! GEMMAT is the value of transa that is used in the first gemm call
         ! GEMMBT is the value of transb that is used in the first gemm call
         CHARACTER         TRMMS, TRMMT, TRMM1U, TRMM2U, GEMMAT, GEMMBT
         ! Parameters
         DOUBLE PRECISION ONE, NEG_ONE, ZERO
         PARAMETER(ONE=1.0D+0, ZERO = 0.0, NEG_ONE=-1.0D+0)
         ! We change the algorithm depending on the values of DIRECT and STOREV
         ! If DIRECT='F', then  H = H(1)...H(k)
         ! Otherwise,           H = H(k)...H(1)
         ! If STOREV='C', then  the reflectors are stored as columns in V and T
         !      will be upper triangular.
         ! Otherwise, the reflectors are stored as rows in V and T is lower
         !      triangular
         ! If 
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

         ! Beginning of executable statements
         ! Early exit if possible
         IF(N.EQ.0) THEN
            RETURN
         END IF
         ! Base case
         IF(N.EQ.1) THEN
            T(1,1) = TAU(1)
            RETURN
         END IF
         ! Compute necessary flags and indices before doing computations
         ! TODO: Consider having a subroutine that computes these flags once and
         ! passes them into a 'recursive' subroutine
         ! K = floor(N/2)
         K = N / 2
         ! Determine if we are going in the forward or backward direction
         DIRF = LSAME(DIRECT,'F')
         ! Determine if the reflectors are stored as column or row vectors
         COLV = LSAME(STOREV,'C')
         ! Determine if we copy over V1^\top
         ! This happens when (STOREV='C' and DIRECT='F') or (STOREV='R' and DIRECT='B')
         COPYTR = (DIRF.AND.COLV).OR.((.NOT.COLV).AND.(.NOT.DIRF))
         ! Determine if we will multiply by the transpose of the triangular
         ! matrix. This will happen when .NOT.COPYTR
         TRMMT = 'T'
         IF(COPYTR) THEN
            TRMMT = 'N'
         END IF
         ! Determine transpose flags for call to GEMM
         GEMMAT = 'N'
         GEMMBT = 'T'
         ! Determine if the first triangular call will be upper or lower
         TRMM1U = 'U'
         ! Determine if we are storing T as upper or lower
         TRMM2U = 'L'
         ! Determine if we want to put T1 or T2 on the left when we finish
         ! computing T3
         TLEFTI = K+1
         TLEFTJ = K+1
         TRIGHTI= 1
         TRIGHTJ= 1
         IF(COLV) THEN
            GEMMAT = 'T'
            GEMMBT = 'N'
            TRMM1U = 'L'
            TRMM2U = 'U'
            TLEFTI = 1
            TLEFTJ = 1
            TRIGHTI= K+1
            TRIGHTJ= K+1
         END IF
         ! Determine the indices for the 'A' and 'B' matrices in the call to
         ! GEMM
         AI = N+1
         AJ = K+1
         BI = N+1
         BJ = 1
         ! Determine what character flag should be used for the call to trmm for
         ! the side of the triangular matrix
         TRMMS = 'L'
         IF(DIRF) THEN
            AI = N+1
            AJ = 1
            BI = N+1
            BJ = K+1
            TRMMS = 'R'
         END IF
         ! If V is stored as rows, then V1 starts at row 1 and column K+1
         V1I = 1
         V1J = K+1
         ! If V is stored as columns, then V1 starts at row K+1 and column 1
         IF(COLV) THEN
            V1I = K+1
            V1J = 1
         END IF
         ! Determine where T_3 lives
         ! If we are going backward, then T_3 is in the lower triangular part
         T3I = K+1
         T3J = 1
         IF(DIRF) THEN
            ! If we are going forward, then T_3 is in the upper triangular part
            T3I = 1
            T3J = K+1
         END IF
         ! Compute T_1
         CALL MY_DLARFT_REC(DIRECT, STOREV, M, K, V, LDV, TAU, T, LDT)

         ! Compute T_2
         CALL MY_DLARFT_REC(DIRECT, STOREV, M-K, N-K, V(K+1,K+1),
     $         LDV, TAU(K+1), T(K+1,K+1), LDT)

         ! Compute T_3 
         ! First, we compute V1^\top * V2 or V1 * V2\top
         ! Copy in V1 or V1^\top
         IF(COPYTR) THEN
            ! Copying in V1^\top manually as we don't have a routine to copy transposes
            ! Below may not properly work for when DIRECT='B'
            DO I = 1, K
               DO J = 1, N-K
                  T(T3I+I-1,T3J+J-1) = V(V1I+J-1, V1J+I-1)
               END DO
            END DO
         ELSE
            ! Copy in V1
            CALL DLACPY('All', N-K, K, V(V1I,V1J), LDV, T(T3I,T3J), LDT)
         END IF
         ! T_3 = V_{2,1}^\top * V_{2,2}
         CALL DTRMM(TRMMS, TRMM1U, TRMMT, 'Unit', 
     $         K, N - K, ONE, V(K+1, K+1), LDV, T(T3I, T3J), LDT)

         IF(M.GT.N) THEN
         ! T_3 = T_3 + V_{3,1}^\topV_{3,2}
            CALL DGEMM(GEMMAT, GEMMBT, K, N-K, M-N, ONE,
     $            V(N+1, 1), LDV, V(N+1,K+1), LDV, ONE, 
     $            T(T3I, T3J), LDT)
         END IF

         ! At this point, we have that T_3 = V_1^\top *V_2
         ! All that is left is to pre and post multiply by -T_1 and T_2
         ! respectively.

         ! T_3 = -T_1*T_3
         CALL DTRMM('Left', TRMM2U, 'No transpose', 'Non-unit',
     $         K, N - K, NEG_ONE, T(TLEFTI,TLEFTJ), LDT, T(T3I, T3J), 
     $         LDT)
         ! T_3 = T_3*T_2
         CALL DTRMM('Right', TRMM2U, 'No transpose', 'Non-unit',
     $         K, N - K, ONE, T(TRIGHTI,TRIGHTJ), LDT, T(T3I,T3J), LDT)
         
         ! Now, we have T in the correct form
      END SUBROUTINE
