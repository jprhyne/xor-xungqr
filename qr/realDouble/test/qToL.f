*     This routine copies L into the matrix X where Q is the result of 
*     calling dgeqlf
      SUBROUTINE QTOL(M, N, Q, X)
         ! Arguments
         INTEGER M, N
         DOUBLE PRECISION Q(M,*), X(N,*)

         ! Local variables
         INTEGER I, J

         DO I = 1, N
            DO J = 1, N
               X(I,J) = Q(M-N+I, J)
            END DO
         END DO
      END SUBROUTINE
