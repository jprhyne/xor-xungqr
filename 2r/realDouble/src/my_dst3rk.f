c Recursively computes the rank K symmetric update of the form
c  C = V**T * V 
c  where C is a symmetric matrix, and V is unit lower triangular 
c   This is an im place version of the existing out of place version.
c Note:
c Future support for computing V*V**T or for V being non-unit or V being upper
c  triangular is possible, but not currently planned unless comes up as a
c  potential use case in our development.
c Arguments
c  M (in)      - Number of rows in V and M\geq N
c  N (in)      - Number of rows of C as well as the number of columns of V
c  V (in/out)  - M by M matrix that has assumed ones on the diagonal and is
c                 assumed to contain only 0's above the diagonal.
c              On output: the upper and diagonal part of V will contain the
c                 product.
c  LDV (in)    - leading dimension of V
c     Cost: (2m^3 + 3m^2 -2m)/6
      RECURSIVE SUBROUTINE MY_DST3RK(M, V, LDV)
         ! Scalar arguments
         INTEGER           M, LDV
         ! Array arguments
         DOUBLE PRECISION  V(LDV,*)
         ! Local variables
         INTEGER           K, INFO, I, J
         ! Parameters
         DOUBLE PRECISION  ONE, ZERO
         PARAMETER(ONE=1.0D+0, ZERO=0.0D+0)
         ! External functions
         EXTERNAL DSYRK, DTRMM
         ! We break down V as follows (where K = N/2) [This is because we assume
         ! V is square]
         ! V = |---------|
         !     |V_11 0   |
         !     |V_21 V_22|
         !     |---------|
         ! V_11 \in \R^{K\times K}
         ! V_21 \in \R^{M-K\times K}
         ! V_22 \in \R^{M-K\times M-K}
         ! Writing out the multiplication using this break down, we get
         ! V**T * V = |---------| ** T |---------|
         !            |V_11 0   |      |V_11 0   |
         !            |V_21 V_22|      |V_21 V_22|
         !            |---------|      |---------|
         !          = |---------------||---------|
         !            |V_11**T V_21**T||V_11 0   |
         !            |0       V_22**T||V_21 V_22|
         !            |---------------||---------|
         !          = |----------------| + |------------------------------|
         !            |V_11**T * V_11 0|   |V_21**T * V_21 V_21**T * V_22 |
         !            |0              0|   |V_22**T * V_21 V_22**T * V_22 |
         !            |----------------|   |------------------------------|
         ! However since the result is a symmetric matrix and we are doing our
         ! operations in place, we omit computing (V**T * V)_21. This leaves us
         ! with only needing to compute the following pieces
         ! V_11 = V_11**T * V_11 + V_21**T * V_21
         ! V_12 = V_21**T * V_22          [TRMM after a copy]
         ! V_22 = V_22**T * V_22          [This routine]
         ! V_11 can be broken down further to get
         ! V_11 = V_11**T * V_11          [This routine]
         ! V_11 = V_11 + V_21**T * V_21   [SYRK]
         ! 
         ! Beginning of executable statements
         ! Quick exit if possible
         IF( M.EQ.0 ) THEN
            RETURN
         END IF
         ! Checking for terminating case
         IF( M.EQ.1 ) THEN
            V(1,1) = ONE
            RETURN
         END IF
         ! Now, we are in the recursive case
         ! Compute K
         K = M / 2
         ! Compute V_11
         CALL MY_DST3RK(K, V(1,1), LDV)
         CALL DSYRK('U', 'T', K, M-K, ONE, V(K+1,1), LDV, ONE, 
     $               V(1,1), LDV)
         ! Compute C_12
         ! Copy over V_21**T. This may be able to be moved to be the user's
         ! responsibility before calling this
         ! Since we are copying a transpose of a matrix and no existing routines
         ! exist to do this (that I've found) we do this manually
         DO I = K+1, M
            DO J = 1, K
               V(J,I) = V(I,J)
            END DO
         END DO
         CALL DTRMM('R', 'L', 'N', 'U', K, M-K, ONE, V(K+1, K+1), LDV, 
     $               V(1, K+1), LDV)
c         CALL DTRMMOOP('R', 'U', K, N-K, V(K+1, 1), LDV, V(K+1,K+1),
c     $               LDV, BETA, C(1, K+1), LDC, INFO)
         ! Compute C_22
         CALL MY_DST3RK(M-K, V(K+1, K+1), LDV)

      END SUBROUTINE
