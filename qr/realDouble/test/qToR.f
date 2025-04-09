*     This routine copies R into the matrix X where Q is the result of 
*     calling dgerqf
      SUBROUTINE QTOR(M, N, Q, X)
         ! Arguments
         INTEGER M, N
         DOUBLE PRECISION Q(M,*), X(M,*)

         ! Local variables
         INTEGER I, J

         DO I = 1, M
            DO J = 1, M
               X(I,J) = Q(I, N-M+J)
            END DO
         END DO
      END SUBROUTINE
