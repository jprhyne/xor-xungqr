*     Cost: (2mn**2 + n**2 - n)/2
      SUBROUTINE MY_ZORGKR(M, N, Q, LDQ)
         ! Arguments
         INTEGER           M, N, LDQ

         ! Array arguments
         COMPLEX*16  Q(LDQ,*)

         ! Scalar variables
         INTEGER           I, J

         ! External subroutines
         EXTERNAL ZTRMM, ZTVT, ZLUMM

         ! Parameters
         COMPLEX*16        NEG_ONE
         PARAMETER(NEG_ONE=(-1.0D+0,0.0D+0))
         

         CALL ZTVT(N, Q, LDQ)

         ! Now, we have computed T=TV_1^\top

         ! Compute Q = -VT
         ! first n rows
         ! Q_1 = -V_1*T ( LUMM lower-upper matrix mult )
         ! Q_2 = -V_2*T (TRMM)
         
         ! Computing Q_2 only when necessary
         IF (M.GT.N) THEN
            CALL ZTRMM('Right', 'Upper', 'No-transpose', 'Non-unit',
     $         M-N, N, NEG_ONE, Q, LDQ, Q(N+1,1), LDQ)
         END IF

         ! Compute -Q_1
         CALL ZLUMM(N, NEG_ONE, Q, LDQ)

         ! Compute "I" - Q
         J = MIN(M,N)
         DO I = 1, J
            Q(I,I) = Q(I,I) + 1.0D+0
         END DO

         ! Now, we should have Q that satisfies the conditions of org2r. IE Q
         ! has orthonormal columns and is a Q associated with A in the QR
         ! decomposition

      END SUBROUTINE
