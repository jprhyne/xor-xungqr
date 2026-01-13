      PROGRAM COUNT_LUMM
         INTEGER           N, LDA
         CHARACTER         SIDEL, DIAGL, DIAGU
         DOUBLE PRECISION  ALPHA, COUNT, MY_VAL, N_DBLE, MY_COUNT

         INTEGER           I

         INTEGER, DIMENSION(4) :: NVEC
         DOUBLE PRECISION, ALLOCATABLE :: A(:,:)
         NVEC(1) = 3
         NVEC(2) = 5
         NVEC(3) = 12
         NVEC(4) = 123

         SIDEL = 'L'
         DIAGL = 'N'
         DIAGU = 'U'
         ALPHA = 1.0D+0
         LDA = NVEC(4)

         ALLOCATE(A(LDA,LDA))

         CALL RANDOM_NUMBER(A)

         DO I = 1, 4
            N = NVEC(I)
            COUNT = 0.0D+0
            CALL DLUMM_CNT(SIDEL, DIAGL, DIAGU, N, ALPHA, A, LDA, COUNT)
            N_DBLE = DBLE(N)
            MY_VAL = MY_COUNT(N_DBLE)
            WRITE(*,*) "N = ", N, "COUNT = ", COUNT, "MY_VAL = ", MY_VAL
         END DO

         DEALLOCATE(A)
      END PROGRAM

      DOUBLE PRECISION FUNCTION MY_COUNT(N)
         DOUBLE PRECISION :: N

         MY_COUNT = N*N*N - N
         MY_COUNT = 2.0 * MY_COUNT / 3.0 + N
      END FUNCTION
