      PROGRAM TEST_ZLARFT_UT
         ! This program tests that zlarft_ut produces a T that satisfies the
         ! same constraints as zlarft_lvl2
         ! Local variables
         INTEGER           M, N, LWORK, INFO
         INTEGER :: SEEDVEC(4)
         DOUBLE PRECISION  NORM_ORTH, NORM_REP, NORM_A
         ! Local arrays
         COMPLEX*16, ALLOCATABLE :: A(:,:), Q(:,:), T(:,:), R(:,:),
     $      V(:,:), WORK(:,:), TAU(:)
         ! External Subroutines
         EXTERNAL          ZHTRK, ZLACPY, ZLASET
         ! External Functions
         DOUBLE PRECISION  ZLANGE
         EXTERNAL          ZLANGE
         ! Parameters
         COMPLEX*16        ONE, ZERO, NEG_ONE
         PARAMETER(ONE=(1.0D+0,0.0D+0), ZERO=(0.0D+0,0.0D+0),
     $      NEG_ONE=(-1.0D+0,0.0D+0))
         ! Intrinsic functions
         INTRINSIC         SQRT,DBLE
         ! Query the user for the problem dimensions
         WRITE(*,*) "Provide an M and N for A (M.GE.N)"
         READ(*,*) M, N
         ! We will not do error checking, just hope the user provided the
         ! correct relation of M and N. Instantiate the seed array. 
         SEEDVEC(1) = 457
         SEEDVEC(2) = 805
         SEEDVEC(3) = 646
         SEEDVEC(4) = 789

         ! Allocate our arrays
         ALLOCATE(A(M,N))
         ALLOCATE(V(M,N))
         ALLOCATE(Q(M,N))
         ALLOCATE(T(N,N))
         ALLOCATE(R(N,N))
         ALLOCATE(WORK(M,M)) ! This is an over-allocation... too bad!
         ALLOCATE(TAU(M))
         ! Generate A as a random matrix
         CALL ZLARNV(4, SEEDVEC, M*N, A)
         ! Compute the norm of A for later use
         NORM_A = ZLANGE('F', M, N, A, M, WORK)
         ! Copy A into V
         CALL ZLACPY('All', M, N, A, M, V, M)
************************************************************************
*        QR decomposition                                              *
************************************************************************
         LWORK = INT(M*M)
         ! Factorize V
         CALL ZGEQRF(M, N, V, M, TAU, WORK, LWORK, INFO)
         ! Copy V into Q
         CALL ZLACPY('All', M, N, V, M, Q, M)
         ! Copy R into R
         CALL ZLACPY('Upper', N, N, V, M, R, N)
         ! Form T using the level2 implementation
         CALL ZLARFT_LVL2('F', 'C', M, N, Q, M, TAU, T, N)
         ! Copy T into the upper triangular part of Q
         CALL ZLACPY('Upper', N, N, T, N, Q, M)
         ! Call dorgkr to form Q
         CALL ZUNGKR(M, N, Q, M)
         ! Initialize WORK as the identity matrix
         CALL ZLASET('All', M, M, ZERO, ONE, WORK, M)
         ! Now we test that Q is orthogonal
         CALL ZHERK('U', 'C', N, M, ONE, Q, M, NEG_ONE, WORK, M)
         ! Compute the frobenius norm of Q**H*Q - I
         NORM_ORTH = ZLANGE('F', N, N, WORK, M, WORK)
         NORM_ORTH = NORM_ORTH / SQRT(DBLE(N))
         ! Copy A into Work
         CALL ZLACPY('All', M, N, A, M, WORK, M)
         ! Compute WORK = Q*R - WORK
         CALL ZTRMMOOP('R','U','N','N','N', M, N, ONE, R, N, Q, M, 
     $      NEG_ONE, WORK, M)
         ! Compute the norm of QR-A
         NORM_REP = ZLANGE('F', M, N, WORK, M, WORK)
         NORM_REP = NORM_REP / NORM_A
         WRITE(*,*) "ZLARFT_LVL2"
         WRITE(*,*) "||Q**H * Q - I||_F/||I||_F = ", NORM_ORTH
         WRITE(*,*) "||Q*R - A||_F / ||A||_F    = ", NORM_REP
         ! Copy V into Q
         CALL ZLACPY('All', M, N, V, M, Q, M)
         ! Fom T using the UT based implementation
         CALL ZLARFT_UT('F', 'C', '1', M, N, Q, M, TAU, T, N)
         ! Copy T into the upper triangular part of Q
         CALL ZLACPY('Upper', N, N, T, N, Q, M)
         ! Call dorgkr to form Q
         CALL ZUNGKR(M, N, Q, M)
         ! Initialize WORK as the identity matrix
         CALL ZLASET('All', M, M, ZERO, ONE, WORK, M)
         ! Now we test that Q is orthogonal
         CALL ZHERK('U', 'C', N, M, ONE, Q, M, NEG_ONE, WORK, M)
         ! Compute the frobenius norm of Q**H*Q - I
         NORM_ORTH = ZLANGE('F', N, N, WORK, M, WORK)
         NORM_ORTH = NORM_ORTH / SQRT(DBLE(N))
         ! Copy A into Work
         CALL ZLACPY('All', M, N, A, M, WORK, M)
         ! Compute WORK = Q*R - WORK
         CALL ZTRMMOOP('R','U','N','N','N', M, N, ONE, R, N, Q, M, 
     $      NEG_ONE, WORK, M)
         ! Compute the norm of QR-A
         NORM_REP = ZLANGE('F', M, N, WORK, M, WORK)
         NORM_REP = NORM_REP / NORM_A
         WRITE(*,*) "ZLARFT_UT"
         WRITE(*,*) "||Q**H * Q - I||_F/||I||_F = ", NORM_ORTH
         WRITE(*,*) "||Q*R - A||_F / ||A||_F    = ", NORM_REP
         ! Free our memory
10       DEALLOCATE(A)
         DEALLOCATE(V)
         DEALLOCATE(Q)
         DEALLOCATE(T)
         DEALLOCATE(R)
         DEALLOCATE(WORK)
         DEALLOCATE(TAU)
      END PROGRAM
