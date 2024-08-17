c This file tests the functionality of my_dst3rk and my_dst3rk_oop
c The main difference is the former is in place while the latter is out of place

      SUBROUTINE TEST_DST3RK(M, N) 
         ! Scalar arguments
         INTEGER M, N
         ! Local variables
         ! Scalar variables
         INTEGER I, J
         DOUBLE PRECISION NORM_F, TMP
         ! Array variables
         DOUBLE PRECISION, ALLOCATABLE :: V(:,:), Vs(:,:), C(:,:)
         ! Parameters
         DOUBLE PRECISION ONE, ZERO
         PARAMETER(ONE=1.0D+0, ZERO=0.0D+0)
         ! First we test the in place operation as it is the easiest to do
         ! Allocate memory
         ALLOCATE(V(M,M))
         ALLOCATE(Vs(M,M))
         ALLOCATE(C(M,M))
         ! Generate V as a random matrix
         CALL RANDOM_NUMBER(V)
         ! Store V inside Vs for checking that we don't touch the lower part
         CALL DLACPY('Lower', M, M, V, M, Vs, M)
         ! Call my in place version
         CALL MY_DST3RK(M, V, M)
         ! Ensure the below diagonal was not touched
         DO I = 2, M
            DO J = 1, I-1
               IF( V(I,J).NE.Vs(I,J) ) THEN
                  WRITE(*,*) "Inconsistency in V at I = ", I, 
     $                           " and J = ", J
               END IF
            END DO
         END DO
         ! For the current implementation of DSYRK, we must manyally set the
         ! upper part of Vs to be 0 and the diagonal to be 1
         DO I = 1, M
            DO J = I, M
               Vs(I,J) = ZERO
            END DO
            Vs(I,I) = ONE
         END DO
         ! Call existing functionality that is out of place. 
         CALL DSYRK('Upper', 'Transpose', M, M, ONE, Vs, M, ZERO, C, M)

         NORM_F = 0.0
         DO I = 1, M
            DO J = I, M
               TMP = C(I,J) - V(I,J)
               NORM_F = NORM_F + TMP * TMP
            END DO
         END DO
         ! Print out the frobenius norm of the differnce in C and V
         WRITE(*,*) "||C - V||_F = ", NORM_F

         DEALLOCATE(V)
         DEALLOCATE(Vs)
         DEALLOCATE(C)

      END SUBROUTINE 
