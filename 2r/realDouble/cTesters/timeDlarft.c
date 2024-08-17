#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
double computeRecPerf(double time, double m, double n)
{
    double rec = ((double) n * (double) n - 1.0) * (2.0 * (double) m + (double) n);
    rec /= 6.0;
    return rec / (time*1.0e+9);
}
double computeUTPerf(double time, double m, double n) 
{
    //  (2m^3 + 3m^2 + 6mn^2 + 6mn - 5m + 2n^3 -12n^2 + 10n)/6 
    double ut = 2.0*m*m*m + 3.0*m*m + 6.0*m*n*n + 6.0*m*n - 5.0*m + 2.0*n*n*n - 12.0*n*n + 10.0*n;
    ut /= 6.0;
    return ut / (time*1.0e+9);
}
// Compute the orthogonality norm
// ||Q**T * Q - I||
double computeOrthNorm(int m, int n, double *Q)
{
    double *work = (double *) malloc(n * n * sizeof(double));
    char aChar = 'a';
    char uChar = 'u';
    char tChar = 't';
    char fChar = 'f';

    double one    =  1.0;
    double zero   =  0.0;
    double negOne = -1.0;
    int dummy = 0;
    
    // Set work to be I
    dlaset_(&aChar, &n, &n, &zero, &one, work, &n, dummy);
    // Compute work = Q**T * Q - I
    dsyrk_(&uChar, &tChar, &n, &m, &one, Q, &m, &negOne, work, &n);
    // Compute the norm of work
    double ret = 0.0;
    for (int i = 0; i < n * n; i++)
        ret += work[i] * work[i];
    ret = sqrt(ret);
    return ret;
}
// Compute the representation norm
// ||As - QR|| / ||As||
double computeRepresNorm(int m, int n, double *Q, double *R, int ldr, double *As, int lda, double normA)
{
    char aChar = 'a';
    char rChar = 'r';
    char uChar = 'u';
    char nChar = 'n';
    char fChar = 'f';
    int dummy = 0;
    double one = 1.0;
    double *work = (double *)malloc(m * n * sizeof(double));
    // Copy Q into work
    dlacpy_(&aChar, &m, &n, Q, &m, work, &m, dummy);

    // Compute work = work * R
    dtrmm_(&rChar, &uChar, &nChar, &nChar, &m, &n, &one, R, &ldr, work, &m, dummy, dummy, dummy, dummy);
    
    // Compute work = work - A
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            work[i + j * m] -= As[i + j * lda];
    double ref = 0.0;
    for (int i = 0; i < m * n; i++) 
        ref += work[i]*work[i];
    ref = sqrt(ref);
    return ref / normA;
}
// In order to ensure accuracy, we will throw the computed T factor into 
// my_dorgkr, and checking the relative error in computing Q and R.
// This does implicitly rely on my_dorgkr working properly. So if the call
// to either reference or optimized dlartf fails, first ensure that routine
// is behaving properly.
int main(int argc, char *argv[]) {
    int info, lda, ldq, m, n, lwork, nb, i, j;
    int workQuery = -1;
    // double variables
    double *A, *As, *V, *Vs, *T, *work, *tau=NULL;
    double normA;
    double perform_refL, elapsed_refL;
    double norm_orth_ref, norm_repres_ref;
    double norm_orth_opt, norm_repres_opt;
    double norm_orth_rec, norm_repres_rec;
    double norm_orth_ut, norm_repres_ut;
    double zero = 0;
    double one = 1;
    double negOne = -1;
    // struct for help with timing
    struct timeval tp;
    // character variables
    char aChar = 'A';
    char fChar = 'F';
    char lChar = 'L';
    char nChar = 'N';
    char rChar = 'R';
    char tChar = 'T';
    char uChar = 'U';
    char cChar = 'C';

    // Dummy value that is used when calling fortran
    // functions that take characters from C
    size_t dummy = 0;

    bool timesOnly = false;

    m = 30;
    n = -1;
    lda = -1;
    ldq = -1;

    for(i = 1; i < argc; ++i){
        if( strcmp( *(argv + i), "-ldq") == 0) {
            ldq  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-lda") == 0) {
            lda  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-m") == 0) {
            m  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-n") == 0) {
            n  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-t") == 0) {
            timesOnly = true;
        }
    }
    //ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
    char *dorgqr = "dorgqr";
    char *space = "";
    int intOne = 1;
    int intNegOne = -1;
    // Determines blocking parameter for a large value of m,n, and k. This is to imitate a single panel
    // By calling the ilaenv we are compiling against,
    // we can ensure that even though we don't know exactly what ilaenv is doing, we are using the same blocking parameter
    // We only call this routine if the user does not provide a value for n
    if (n == -1)
        n = ilaenv_( &intOne, dorgqr, space, &m, &m, &m, &intNegOne, dummy, dummy);
    // In the event n was too big, we force n to be m. This is because we do not do proper
    // error checking for sake of speed as we are assuming we only working on tall and skinny matrices
    if (n > m)
        n = m;
    if( lda < 0 ) lda = m;
    // allocate memory for the matrices and vectors that we need to use
    A = (double *) malloc( lda * n * sizeof(double));
    As = (double *) malloc( lda * n * sizeof(double));
    V = (double *) calloc( m * n, sizeof(double));
    Vs= (double *) malloc( m * n * sizeof(double));
    T = (double *) malloc( n * n * sizeof(double));

    // Print to the user what we are doing along with any arguments that are used
    if(!timesOnly) {
        printf("dgeqrf dlarft | m = %4d, n = %d, lda = %4d, ldq = %4d\n", m, n, lda, ldq);
    }

    // Generate A as a random matrix
    for (i = 0; i < lda * n; i++)
        A[i] = (double) rand() / (double) (RAND_MAX) - 0.5e+00;

    normA = 0.0;
    for (i = 0; i < lda * n; i++) 
        normA += A[i] * A[i];
    normA = sqrt(normA);
    // Create the work array to do workspace queries
    work = (double *) malloc(sizeof(double));
    // allocate the tau vector
    tau = (double *) malloc(n * sizeof(double));
    // Determine how much workspace is needed for our operations
    dgeqrf_(&m, &n, A, &lda, tau, work, &workQuery, &info );
    lwork = work[0];

    // reallocate work to be of the right size
    work = (double *) realloc(work, lwork * sizeof(double));

    // Store a copy of A to use later
    dlacpy_(&aChar, &m, &n, A, &lda, As, &lda, dummy);

    // Call dgeqrf
    dgeqrf_(&m, &n, A, &lda, tau, work, &lwork, &info);
    // At this point, we now have that A = [R \\ V]
    // Copy the 'lower triangular' part of A into V
    dlacpy_(&lChar, &m, &n, A, &lda, V, &m, dummy);
    // Set the diagonal of V to be exactly 1.
    for( i = 0; i < n; i++)
        V[i+i*m] = one;
    // Now, V is the exact matrix of our householder vectors
    // Store a copy of V into Vs
    dlacpy_(&aChar, &m, &n, V, &m, Vs, &m, dummy);


    // Start the timer
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    // Call reference dlarft using optimized blas as the backend
    dlarft_ref_(&fChar, &cChar, &m, &n, V, &m, tau, T, &n);
    // grab the execution time
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    // Store this value 
    double refTime = elapsed_refL;
    // Compute the number of flops
    double refFlop = computeRecPerf(refTime, (double) m, (double) n);

    // testing if T is valid
    // Copy T into the upper triangular part of V
    dlacpy_(&uChar, &m, &n, T, &n, V, &m, dummy);
    // Call my_dorgkr with the above V matrix
    my_dorgkr_(&m, &n, V, &m);
    // Compute \|Q**T *Q - I\|/\|As\|
    norm_orth_ref = computeOrthNorm(m, n, V);
    // Compute ||A - Q*R||_F / ||A||_F
    norm_repres_ref = computeRepresNorm(m, n, V, A, lda, As, lda, normA);

    // copy Vs back into V in the event V was changed
    dlacpy_(&aChar, &m, &n, Vs, &m, V, &m, dummy);

    
    // Start the timer
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    // Call reference dlarft using optimized blas as the backend
    dlarft_(&fChar, &cChar, &m, &n, V, &m, tau, T, &n, dummy, dummy);
    // grab the execution time
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    // Store this value 
    double optTime = elapsed_refL;
    // Compute the number of flops
    double optFlop = computeRecPerf(optTime, (double) m, (double) n);

    // testing if T is valid
    // Copy T into the upper triangular part of V
    dlacpy_(&uChar, &m, &n, T, &n, V, &m, dummy);
    // Call my_dorgkr with the above V matrix
    my_dorgkr_(&m, &n, V, &m);
    // Compute \|Q**T *Q - I\|/\|As\|
    norm_orth_opt = computeOrthNorm(m, n, V);
    // Compute ||A - Q*R||_F / ||A||_F
    norm_repres_opt = computeRepresNorm(m, n, V, A, lda, As, lda, normA);

    // copy Vs back into V in the event V was changed
    dlacpy_(&aChar, &m, &n, Vs, &m, V, &m, dummy);


    // Start the timer
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    // Call reference dlarft using optimized blas as the backend
    my_dlarft_rec_(&m, &n, V, &m, tau, T, &n);
    // grab the execution time
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    // Store this value 
    double recTime = elapsed_refL;
    // Compute the number of flops
    double recFlop = computeRecPerf(recTime, (double) m, (double) n);

    // testing if T is valid
    // Copy T into the upper triangular part of V
    dlacpy_(&uChar, &m, &n, T, &n, V, &m, dummy);
    // Call my_dorgkr with the above V matrix
    my_dorgkr_(&m, &n, V, &m);
    // Compute \|Q**T *Q - I\|/\|As\|
    norm_orth_rec = computeOrthNorm(m, n, V);
    // Compute ||A - Q*R||_F / ||A||_F
    norm_repres_rec = computeRepresNorm(m, n, V, A, lda, As, lda, normA);

    // copy Vs back into V in the event V was changed
    dlacpy_(&aChar, &m, &n, Vs, &m, V, &m, dummy);


    // Start the timer
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    // Call reference dlarft using optimized blas as the backend
    my_dlarft_ut_(&m, &n, V, &m, tau, T, &n);
    // grab the execution time
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    // Store this value 
    double utTime = elapsed_refL;
    // Compute the number of flops
    double utFlop = computeUTPerf(utTime, (double) m, (double) n);

    // testing if T is valid
    // Copy T into the upper triangular part of V
    dlacpy_(&uChar, &m, &n, T, &n, V, &m, dummy);
    // Call my_dorgkr with the above V matrix
    my_dorgkr_(&m, &n, V, &m);
    // Compute \|Q**T *Q - I\|/\|As\|
    norm_orth_ut = computeOrthNorm(m, n, V);
    // Compute ||A - Q*R||_F / ||A||_F
    norm_repres_ut = computeRepresNorm(m, n, V, A, lda, As, lda, normA);

    // Now, we print out the testing information
    printf("reference DLARFT\ntime: %10.10e\nperf: %10.10e\north: %10.10e\nrepres: %10.10e\n", refTime, refFlop, norm_orth_ref, norm_repres_ref);
    printf("optimized DLARFT\ntime: %10.10e\nperf: %10.10e\north: %10.10e\nrepres: %10.10e\n", optTime, optFlop, norm_orth_opt, norm_repres_opt);
    printf("MY_DLARFT_REC\ntime: %10.10e\nperf: %10.10e\north: %10.10e\nrepres: %10.10e\n", recTime, recFlop, norm_orth_rec, norm_repres_rec);
    printf("MY_DLARFT_UT\ntime: %10.10e\nperf: %10.10e\north: %10.10e\nrepres: %10.10e\n", utTime, utFlop, norm_orth_ut, norm_repres_ut);
}
