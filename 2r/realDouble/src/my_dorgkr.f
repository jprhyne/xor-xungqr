*     Cost: (2mn**2 + n**2 - n)/2
      SUBROUTINE MY_DORGKR(M, N, Q, LDQ)
         ! Arguments
         INTEGER           M, N, LDQ

         ! Array arguments
         DOUBLE PRECISION  Q(LDQ,*)

         ! Scalar variables
         INTEGER           I, J

         DOUBLE PRECISION, ALLOCATABLE :: V(:,:)
         ! External subroutines
         EXTERNAL DTRMM, DTVT, LUMM

         ! Parameters
         DOUBLE PRECISION NEG_ONE
         PARAMETER(NEG_ONE=-1.0D+0)
         

         CALL DTVT(N, Q, LDQ)

         ! Now, we have computed T=TV_1^\top

         ! Compute Q = -VT
         ! first n rows
         ! Q_1 = -V_1*T ( LUMM lower-upper matrix mult )
         ! Q_2 = -V_2*T (TRMM)
         
         ! Computing Q_2 only when necessary
         IF (M.GT.N) THEN
            CALL DTRMM('Right', 'Upper', 'No-transpose', 'Non-unit',
     $         M-N, N, NEG_ONE, Q, LDQ, Q(N+1,1), LDQ)
         END IF

         ! Compute Q_1
         CALL LUMM(N, NEG_ONE, Q, LDQ)

         ! Compute "I" - Q
         J = MIN(M,N)
         DO I = 1, J
            Q(I,I) = Q(I,I) + 1.0D+0
         END DO

         ! Now, we should have Q that satisfies the conditions of org2r. IE Q
         ! has orthonormal columns and is a Q associated with A in the QR
         ! decomposition

      END SUBROUTINE
