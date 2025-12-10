      PROGRAM TIME_GEMM
         INTEGER M,N,K

         REAL START,FINISH

         DOUBLE PRECISION ONE,PERF,TIME

         DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:), C(:,:)

         INTRINSIC CPU_TIME

         ONE = 1.0D+0


         READ(*,*) M, N, K

         ALLOCATE(A(M,N))
         ALLOCATE(B(N,K))
         ALLOCATE(C(M,N))
         CALL CPU_TIME(START)

         CALL DGEMM('N','N',M,N,K,ONE,A,M,B,N,ONE,C,M)

         CALL CPU_TIME(FINISH)

         TIME = FINISH - START
         PERF = 2.0*DBLE(M)*DBLE(N)*DBLE(K) / (1.0D+9 * TIME)

         WRITE(*,*) PERF

         DEALLOCATE(A)
         DEALLOCATE(B)
         DEALLOCATE(C)

      END PROGRAM
