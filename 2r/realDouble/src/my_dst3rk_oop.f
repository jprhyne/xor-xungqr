c Recursively computes the rank K symmetric update of the form
c  C = V**T * V + \beta C
c  where C is a symmetric matrix, and V is unit lower triangular 
c Note:
c Future support for computing V*V**T or for V being non-unit or V being upper
c  triangular is possible, but not currently planned unless comes up as a
c  potential use case in our development.
c We are also only referencing the upper triangular part of C.
c Arguments
c  M (in)      - Number of rows in V and M\geq N
c  N (in)      - Number of rows of C as well as the number of columns of V
c  V (in)      - M by N matrix that has assumed ones on the diagonal and is
c                 assumed to contain only 0's above the diagonal. Only the below
c                 diagonal is referenced
c  LDV (in)    - leading dimension of V
c  BETA (in)   - Scaling factor associated with the C matrix. If BETA is 0, then
c                 the existing elements of C are not referenced. 
c  C (in/out)  - On input, contains an upper triangular matrix that is assumed
c                 to represent an N by N symmetric matrix. Existing elements are
c                 only referenced when BETA is nonzero. If BETA is 0, then C is
c                 output only
c                On output, C will contain the product ALPHA * V**T * V + BETA * C
c  LDC (in)    - Leading dimension of C
      RECURSIVE SUBROUTINE MY_DST3RK_OOP(M, N, V, LDV, BETA, C, LDC)
         ! Scalar arguments
         INTEGER           M, N, LDV, LDC
         DOUBLE PRECISION  ALPHA, BETA
         ! Array arguments
         DOUBLE PRECISION  V(LDV,*), C(LDC,*)
         ! Local variables
         INTEGER           K, INFO
         ! Parameters
         DOUBLE PRECISION  ONE, ZERO
         PARAMETER(ONE=1.0D+0, ZERO=0.0D+0)
         ! External functions
         EXTERNAL DSYRK, DTRMMOOP
         ! We break down V as follows (where K = N/2)
         ! V = |---------|
         !     |V_11 0   |
         !     |V_21 V_22|
         !     |---------|
         ! V_11 \in \R^{K\times K}
         ! V_21 \in \R^{M-K\times K}
         ! V_22 \in \R^{M-K\times N-K}
         ! We do something similar for C to have
         ! C = |---------|
         !     |C_11 C_12|
         !     |0    C_22|
         !     |---------|
         ! C_11 \in \R^{K\times K}
         ! C_12 \in \R^{K\times N - K}
         ! C_22 \in \R^{N-K\times N-K}
         ! C = V**T * V + BETA * C
         ! C = |---------------| |---------| + BETA * |---------|
         !     |V_11**T V_21**T| |V_11 0   |          |C_11 C_12|
         !     |0       V_22**T| |V_21 V_22|          |0    C_22|
         !     |---------------| |---------|          |---------|
         ! By simplifying this, we can write each component of C as
         ! C_11 = V_11**T * V_11 + V_21**T * V_21 + BETA * C_11
         ! C_12 = V_21**T * V_22 + BETA  * C_12
         ! C_22 = V_22**T * V_22 + BETA  * C_22
         ! Computing C_12 is a call to DTRMMOOP
         ! Computing C_22 is a call to this function
         ! Computing C_11 will be broken down into the following computations
         ! C_11 = V_11**T * V_11 + BETA * C_11 (This function)
         ! C_11 = V_21**T * V_21 + C_11        (DSYRK)
         ! Beginning of executable statements
         ! Checking for terminating case
         IF( N.EQ.1 ) THEN
            IF( M.EQ.1 ) THEN
               C(1,1) = ONE + BETA * C(1,1)
            ELSE
               CALL DSYRK('U', 'T', 1, M-1, ONE, V(2,1), LDV, 
     $                  BETA, C, LDC)
               C(1,1) = ONE + C(1,1)
            END IF
            RETURN
         END IF
         ! Now, we are in the recursive case
         ! Compute K
         K = N / 2
         ! Compute C_11
         CALL MY_DST3RK_OOP(K, K, V(1,1), LDV, BETA, C(1,1), LDC)
         CALL DSYRK('U', 'T', K, M-K, ONE, V(K+1,1), LDV, ONE, 
     $               C(1,1), LDC)
         ! Compute C_12
         CALL DTRMMOOP('R', 'U', K, N-K, V(K+1, 1), LDV, V(K+1,K+1),
     $               LDV, BETA, C(1, K+1), LDC, INFO)
         ! Compute C_22
         CALL MY_DST3RK_OOP(M-K, K, V(K+1, K+1), LDV, BETA, C
     $               (K+1,K+1), LDC)

      END SUBROUTINE
