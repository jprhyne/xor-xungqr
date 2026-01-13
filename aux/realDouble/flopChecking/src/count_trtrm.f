      PROGRAM COUNT_TRTRM
         INTEGER           N, LDT, LDV
         CHARACTER         SIDE, UPLO, TRANSV, DIAGT, DIAGV
         DOUBLE PRECISION  ALPHA, COUNT, MY_VAL, N_DBLE, MY_COUNT

         INTEGER           I

         INTEGER, DIMENSION(4) :: NVEC
         DOUBLE PRECISION, ALLOCATABLE :: T(:,:), V(:,:)
         NVEC(1) = 3
         NVEC(2) = 5
         NVEC(3) = 12
         NVEC(4) = 123

         SIDE = 'R'!ight
         UPLO = 'U'!pper
         TRANSV = 'T'!ranspose
         DIAGT = 'N'!on-unit
         DIAGV = 'U'!nit
         ALPHA = 1.0D+0
         LDT = NVEC(4)
         LDV = NVEC(4)

         ALLOCATE(T(LDT,LDT))
         ALLOCATE(V(LDV,LDV))

         CALL RANDOM_NUMBER(T)
         CALL RANDOM_NUMBER(V)

         DO I = 1, 4
            N = NVEC(I)
            COUNT = 0.0D+0
            CALL DTRTRM_CNT(SIDE, UPLO, TRANSV, DIAGT, DIAGV, N, ALPHA,
     $            T, LDT, V, LDV, COUNT)
            N_DBLE = DBLE(N)
            MY_VAL = MY_COUNT(N_DBLE)
            WRITE(*,*) "N = ", N, "COUNT = ", COUNT, "MY_VAL = ", MY_VAL
         END DO

         DEALLOCATE(T)
         DEALLOCATE(V)
      END PROGRAM

      DOUBLE PRECISION FUNCTION MY_COUNT(N)
         DOUBLE PRECISION :: N

         MY_COUNT = N*N*N + 3.0D+0*N*N - 4.0D+0*N
         MY_COUNT = MY_COUNT / 3.0D+0
         MY_COUNT = MY_COUNT + 2.0D+0*N
      END FUNCTION
