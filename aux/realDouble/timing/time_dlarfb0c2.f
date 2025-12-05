      PROGRAM TIME_LARFB
         ! Local scalars
         INTEGER M, N, K, MK
         REAL START, FINISH
         ! Local Arrays
         DOUBLE PRECISION, ALLOCATABLE :: V(:,:), A(:,:), T(:,:), C(:,:)

         !{TIME,PERF}VEC = [ AOCL (LARFB), REFERENCE (LARFB), LARFB0C2 ]
         DOUBLE PRECISION, DIMENSION(3) :: TIMEVEC, PERFVEC

         INTRINSIC   CPU_TIME

         WRITE(*,*) "Provide M,N,K for LARFB{,0C2}"
         READ(*,*) M, N, K

         MK = M*K

         ALLOCATE(A(MK,MK))
         ALLOCATE(V(MK,MK))
         ALLOCATE(T(K,K))
         ALLOCATE(C(MK,MK))

         CALL RANDOM_NUMBER(A)
         CALL RANDOM_NUMBER(V)
         CALL RANDOM_NUMBER(T)
         CALL RANDOM_NUMBER(C)
      END PROGRAM
