*     Parameters
*     M  :Number of columns in A
*     N  :Number of rows in A
*     K  :Number of reflectors (K\leq N)
*     KI :Column we start at
*     NB :Number of columns in our blocking method
*     A  :Matrix that stores our reflectors on input
*     LDA:Leading dimension of A
*     TAU:Tau vector as a result of geqrf
*
      SUBROUTINE COMPUTE_Q2(M, N, K, KI, NB, A, LDA, TAU)
         ! Arguments
         INTEGER M, N, K, KI, NB, LDA, LDT 
         DOUBLE PRECISION A(LDA, *), TAU(*)
         ! Local variables
         INTEGER I, IB
         DO I = KI + 1, 1, -NB
            IB = NB
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL MY_DLARFT_REC(M-I+1, IB, A(I,I), LDA, TAU(I),
     $         A(I, I), LDA)

*
*           Apply H to A(i:m,i+ib:n) from the left
*
            CALL MY_DLARFB(M - I + 1, N - K, IB, A(I, I), LDA, 
     $         A(I, K + 1), LDA)

         END DO
*        This checks for if K was a perfect multiple of NB
*        so that we only have a special case for the last block when
*        necessary
         IF(I.LT.1) THEN
            IB = I + NB - 1
            I = 1
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL MY_DLARFT_REC(M-I+1, IB, A(I,I), LDA, TAU(I), A(I,I), 
     $         LDA)
*
*           Apply H to A(i:m,i+ib:n) from the left
*
            CALL MY_DLARFB(M - I + 1, N - K, IB, A(I, I), LDA, 
     $         A(I, K + 1), LDA)
         END IF
      END SUBROUTINE
