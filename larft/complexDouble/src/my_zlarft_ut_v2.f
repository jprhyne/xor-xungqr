      SUBROUTINE MY_ZLARFT_UT_V2(M, N, V, LDV, TAU)
         ! Arguments
         ! Scalars
         INTEGER           M, N, LDV
         ! Matrix variables
         COMPLEX*16        V(LDV,*), TAU(N)
         ! Local variables
         INTEGER           I,J,K,INFO
         ! Parameters
         COMPLEX*16        ONE, NEG_ONE
         PARAMETER(ONE = (1.0D+0,0.0D+0))
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
         CALL MY_ZST3RK(N, V, LDV)
         ! Compute V = V_2**T * V_2 + V
         CALL ZSYRK('Upper', 'Transpose', N, M-N, ONE, V(N+1,1), LDV,
     $               ONE, V, LDV)
         ! Sets the diagonal equal to our Tau vector
         DO I = 1, N
            V(I,I) = TAU(I)
         END DO
         ! Replaces V with T^{-1}
         CALL ZTRTRI('Upper', 'Non unit', N, V, LDV, INFO)
      END SUBROUTINE
