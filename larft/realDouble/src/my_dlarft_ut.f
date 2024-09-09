*     Cost:
*     m > n: (6mn^2 + 6mn - 2n^3 - 3n^2 + 2n)/6
*     m = n: (4n^3 - 3n^2 + 2n)/6 
      RECURSIVE SUBROUTINE MY_DLARFT_UT(M, N, V, LDV, TAU)
         ! Arguments
         ! Scalars
         INTEGER           M, N, LDV
         ! Matrix variables
         DOUBLE PRECISION  V(LDV,*), TAU(N)
         ! Local variables
         INTEGER           I,J,K,INFO
         ! Parameters
         DOUBLE PRECISION ONE, NEG_ONE, ZERO, HALF
         PARAMETER(ONE=1.0D+0, HALF=0.5D+0, ZERO = 0.0)
         ! Implementation of the algorithm listed in the following paper
         ! https://www.cs.utexas.edu/users/flame/pubs/p169-joffrain.pdf
         ! Compute T = V^\top V
         ! We do this by breaking up V into the following form
         ! V = |-----|
         !     | V_1 |
         !     | V_2 |
         !     |-----|
         ! where V_1 is unit lower triangular of dimension n by n
         ! V_2 is a rectangular matrix of dimension m-n by n
         ! This means that we can write V**T * V as
         ! V**T = V
         ! |----------------| |-----| = V_1**T * V_1 + V_2**T * V_2
         ! | V_1**T  V_2**T | | V_1 |
         ! |----------------| | V_2 |
         !                    |-----|
         ! Compute this in two parts
         ! V = V_1**T * V_1
         ! Then 
         ! V = V_2**T * V_2 + V
         ! Compute V = V_1**T * V_1
         CALL MY_DST3RK(N, V, LDV)
         ! Compute V = V_2**T * V_2 + V
         CALL DSYRK('Upper', 'Transpose', N, M-N, ONE, V(N+1,1), LDV,
     $               ONE, V, LDV)
         ! Scales the diagonal by 1/2
         CALL DSCAL(N, HALF, V(1,1), LDV + 1)
         ! Replaces V with T^{-1}
         CALL DTRTRI('Upper', 'Non unit', N, V, LDV, INFO)
      END SUBROUTINE
