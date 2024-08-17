      SUBROUTINE TEST_DTVT(N, LDT, LDV)
         ! Arguments
         INTEGER           N, LDT, LDV

         ! Local variables
         INTEGER           I, J, LDQ
         DOUBLE PRECISION  NORM_F, TMP

         ! Local arrays
         DOUBLE PRECISION, ALLOCATABLE :: T(:,:),  V(:,:)
         DOUBLE PRECISION, ALLOCATABLE :: Ts(:,:), Vs(:,:)
         DOUBLE PRECISION, ALLOCATABLE :: A(:,:), Q(:,:)

         ! External subroutines
         EXTERNAL DTVT DLACPY

         ! Parameters
         DOUBLE PRECISION ONE
         PARAMETER(ONE=1.0D+0)

         ! Begin of executable statements
         LDQ = 2*LDT ! for now, we just do this to allow for future testing of the
         ! leading dimension functionality.
         ! Allocate memory
         ALLOCATE(T(LDT,N))
         ALLOCATE(V(LDV,N))
         ALLOCATE(Ts(LDT,N))
         ALLOCATE(Vs(LDV,N))
         ALLOCATE(A(N,N))
         ALLOCATE(Q(LDQ,N))

         ! Generate random T and V
         CALL RANDOM_NUMBER(T)
         CALL RANDOM_NUMBER(V)

         ! Enforce T being upper triangular and V being
         ! unit lower triangular

         DO I = 1, N
            DO J = 1, N
               IF (I.LT.J) THEN
                  V(I,J) = 0
               ELSE IF (I.EQ.J) THEN
                  V(I,J) = 1
               ELSE
                  T(I,J) = 0
               END IF
            END DO
         END DO

         ! Copy T and V into Ts and Vs respectively
         CALL DLACPY('ALL', LDT, N, T, LDT, Ts, LDT)
         CALL DLACPY('ALL', LDV, N, V, LDV, Vs, LDV)

         ! Copy T and V into Q to be of the form that dtvt is expecting. 
         ! IE, layer them on top of each other. Since the diagonal
         ! component will be T, we copy over V first
         CALL DLACPY('Lower', N, N, V, LDV, Q, LDQ)

         ! copy over T. This will overwrite the diagonal with
         ! the elements we want there
         CALL DLACPY('Upper', N, N, T, LDT, Q, LDQ)

         ! Call my subroutine
         CALL DTVT(N, Q, LDQ)

         ! Make sure we don't touch V
         DO I = 2, N 
            DO J = 1, I-1
               IF (V(I,J).NE.Q(I,J)) THEN
                  WRITE(*,*) "Inconsistency at i = ", I, " j = ", J
               END IF
            END DO
         END DO

         ! Call the expected computation (T = TV**T)
         CALL DTRMM('Right', 'Lower', 'Transpose', 'Unit',
     $      N, N, ONE, V, LDV, T, LDT)

         NORM_F = 0
         DO I = 1, N
            DO J = I, N
               TMP = T(I,J) - Q(I,J)
               NORM_F = NORM_F + TMP * TMP
            END DO
         END DO

         WRITE(*,*) "||T_act - T_mine|| = ", NORM_F

         ! Free memory
         DEALLOCATE(T)
         DEALLOCATE(Ts)
         DEALLOCATE(V)
         DEALLOCATE(Vs)
         DEALLOCATE(A)
         DEALLOCATE(Q)

      END SUBROUTINE
