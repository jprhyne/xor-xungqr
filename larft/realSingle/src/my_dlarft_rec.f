c     Cost: m > n: 1/6 * (n^2-1)(2m+n)
c           m = n: 1/2 * (n^3-n)
      RECURSIVE SUBROUTINE MY_DLARFT_REC(DIRECT, STOREV, M, N, V, LDV, 
     $                                   TAU, T, LDT)
         ! Arguments
         ! Scalars
         CHARACTER         DIRECT, STOREV
         ! M is the length of the reflectors
         ! N is the number of reflectors
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
         INTEGER           T3M,T3N,TMP
         LOGICAL           ROWV, DIRB, COPYTR
         ! V2SIDE is the value of side   that is used in the trmm call for v2
         ! V2TRAN is the value of trans  that is used in the trmm call for v2
         ! V2UPLO is the value of uplo   that is used in the trmm call for v2
         ! T3UPLO is the value of uplo   that is used in the trmm call to
         !        finalize computing T3
         ! GEMMAT is the value of transa that is used in the gemm call for the
         !        third component of v1 and v2
         ! GEMMBT is the value of transb that is used in the gemm call for the
         !        third component of v1 and v2
         CHARACTER         V2SIDE, V2TRAN, V2UPLO, T3UPLO, GEMMAT
         CHARACTER         GEMMBT, T1SIDE, T2SIDE
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
         ! If V consists of column vectors,
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
         ! V_{2,2}\in\R^{n-k,n-k} unit lower triangular
         ! V_{3,2}\in\R^{m-n,n-k} rectangular
         ! If V consists of row vectors,
         ! Break V apart into 6 components
         ! V = |-----------------------|
         !     |V_{1,1} V_{1,2} V_{1,3}|
         !     |0       V_{2,2} V_{2,3}|
         !     |-----------------------|
         ! V_{1,1}\in\R^{k,k} unit upper triangular
         ! V_{1,2}\in\R^{k,n-k} rectangular
         ! V_{1,3}\in\R^{k,m-n} rectangular
         ! 
         ! V_{2,2}\in\R^{n-k,n-k} unit upper triangular
         ! V_{2,3}\in\R^{n-k,m-n} rectangular

         ! If DIRECT='F', then
         ! We will construct the T matrix 
         ! T = |---------------| =  |--------|
         !     |T_{1,1} T_{1,2}|    |T_1  T_3|
         !     |0       T_{2,2}|    |0    T_2|
         !     |---------------|    |--------|
         ! If DIRECT='B', then
         ! We will construct the T matrix 
         ! T = |---------------| =  |--------|
         !     |T_{1,1} 0      |    |T_1  0  |
         !     |T_{2,1} T_{2,2}|    |T_3  T_2|
         !     |---------------|    |--------|

         ! T is the triangular factor attained from block reflectors. 
         ! To motivate the structure, consider the product 
         ! Note: This is for the case of DIRECT='F' and STOREV='C'.
         ! The other combinations follow by either swapping the order of 
         ! multiplication for DIRECT='B' and swapping transposition of 
         ! V_1 and V_2 for each term ie all V_1 go to V_1^\top
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
         ! Determine if we are going in the forward or backward direction
         DIRB = LSAME(DIRECT,'B')
         ! Determine if the reflectors are stored as column or row vectors
         ROWV = LSAME(STOREV,'R')
         ! Determine if we copy over V1^\top
         ! This happens when (STOREV='C' and DIRECT='F') or (STOREV='R' and DIRECT='B')
         COPYTR = (DIRB.AND.ROWV).OR.((.NOT.ROWV).AND.(.NOT.DIRB))
         ! K = floor(NumberOfReflectors/2)
         K = N/2
         ! Determine if we will multiply by the transpose of the triangular
         ! matrix. This will happen when .NOT.COPYTR
         V2TRAN = 'T'
         IF(COPYTR) THEN
            V2TRAN = 'N'
         END IF
         ! Compute some helpful indices and flags to make our algorithm work for
         ! all 4 combinations of DIRECT and STOREV
         ! The case of DIRECT = 'F' and STOREV = 'C' will be the values chosen
         ! if none of the 'IF' blocks are hit
         ! 
         ! Indices and sizes first
         ! Determine where the second component of V1 lives
         V1I = K+1
         V1J = 1
         ! Determine where we are going to store T_3
         T3I = 1
         T3J = K+1
         ! Determine the size of T3
         T3M = K
         T3N = N-K
         ! In all cases we need a GEMM call to multiply the third component of
         ! V1 and V2 with each other if the length of the reflectors is greater
         ! than the number of relfectors (IE M > N)
         !
         ! Determine where the 3rd component of V1 lives
         AI  = N+1
         AJ  = 1
         ! Determine where the 3rd component of V2 lives
         BI  = N+1
         BJ  = K+1
         ! Flags next
         ! Determine the side the first component of V2 will be on in our first
         ! multiplication
         V2SIDE = 'R'
         ! Determine if the first component of V2 is upper or lower triangular
         V2UPLO = 'L'
         ! Determine if T3 is stored in the upper or lower component of T (and
         ! consequently if T1 and T2 are upport or lower triangular)
         T3UPLO = 'U'
         ! Determine if the 'A' matrix in our GEMM call is to be transposed or
         ! not
         GEMMAT = 'T'
         ! Determine if the 'B' matrix in our GEMM call is to be transposed or
         ! not
         GEMMBT = 'N'
         ! Determine if T1 will multiply T3 on the right or left
         T1SIDE = 'L'
         ! Determine if T2 will multiply T3 on the right or left
         T2SIDE = 'R'
         ! Change above values if V is stored in rows
         IF(ROWV) THEN
            ! V1 is transposed, so V1{I,J} and A{I,J} swap
            V1I = 1
            V1J = K+1
            AI  = 1
            AJ  = N+1
            ! V2 is transposed, so B{I,J} swap
            BI  = K+1
            BJ  = N+1
            ! Since V1 and V2 are transposed, our transpose flags for 
            ! GEMM are swapped
            GEMMAT = 'N'
            GEMMBT = 'T'
            ! Since V2 is transposed, the first component of V2 is upper
            ! triangular
            V2UPLO = 'U'
         END IF
         ! Change necessary values if we are going 'backwards' (right to left)
         IF(DIRB) THEN
            ! This is a little more tricky. The places of A and B swap so we
            ! need to swap AI with BI and AJ with BJ
            TMP = AI
            AI  = BI
            BI  = TMP
            TMP = AJ
            AJ  = BJ
            BJ  = TMP
            ! T3 is now stored in the lower triangular part of T
            T3I = K+1
            T3J = 1
            ! The size of T3 also changes
            T3M = N-K
            T3N = K
            ! The T1 and T2 are now lower triangular
            T3UPLO = 'L'
            ! The first component of V2 will now be on the left when we
            ! multiply V1 by it
            V2SIDE = 'L'
            ! T1 now multiplies T3 from the right
            T1SIDE = 'R'
            ! T2 now multiplies T3 from the left
            T2SIDE = 'L'
         END IF
         ! Begin actual computation
         ! Compute T_1
         CALL MY_DLARFT_REC(DIRECT, STOREV, M, K, V, LDV, TAU, T, LDT)

         ! Compute T_2
         CALL MY_DLARFT_REC(DIRECT, STOREV, M-K, N-K, V(K+1,K+1),
     $         LDV, TAU(K+1), T(K+1,K+1), LDT)

         ! Compute T_3 = op(V1) * op(V2).
         ! Note: op(.) is either the transpose of the matrix or the matrix
         ! itself depending on the values of DIRECT and STOREV on input. See
         ! above for more details on what exact operations we are doing.
         IF(COPYTR) THEN
            ! Copying in V1^\top manually as we don't have a routine to copy transposes
            ! T3 = (V1)2^\top
            DO I = 1, T3M
               DO J = 1, T3N
                  T(T3I-1+I,T3J-1+J) = V(V1I-1+J, V1J-1+I)
               END DO
            END DO
         ELSE
            ! Copy in (V1)2
            CALL DLACPY('All', T3M, T3N, V(V1I,V1J), LDV, T(T3I,T3J), 
     $                  LDT)
         END IF
         ! Begin computing T_3 = op(V1) * op(V2) 
         ! [op(.) is either transposed or not depending on the values of direct
         ! and storev. See above for more details]
         CALL DTRMM(V2SIDE, V2UPLO, V2TRAN, 'Unit', 
     $         T3M, T3N, ONE, V(K+1, K+1), LDV, T(T3I, T3J), LDT)

         IF(M.GT.N) THEN
            ! If needed, finish the trailing computaiton of op(V1) * op(V2)
            CALL DGEMM(GEMMAT, GEMMBT, T3M, T3N, M-N, ONE,
     $            V(AI,AJ), LDV, V(BI,BJ), LDV, ONE, 
     $            T(T3I, T3J), LDT)
         END IF

         ! At this point, we have that T_3
         ! All that is left is to pre and post multiply by -T_1 and T_2
         ! respectively.

         ! First, T_3 = -T_1*T_3 or -T_3*T_1
         CALL DTRMM(T1SIDE, T3UPLO, 'No transpose', 'Non-unit',
     $         T3M, T3N, NEG_ONE, T(1,1), LDT, T(T3I, T3J), 
     $         LDT)
         ! Next, T_3 = T_3*T_2 or T_2*T_3
         CALL DTRMM(T2SIDE, T3UPLO, 'No transpose', 'Non-unit',
     $         T3M, T3N, ONE, T(K+1,K+1), LDT, T(T3I,T3J), LDT)
         
         ! Now, we have T in the correct form
      END SUBROUTINE
