#ifdef USE_AOCL
    #define SOURCE "AOCL"
#else
    #define SOURCE "REF_LAPACK"
#endif

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

double computePerform_refL(double m, double n, double k, double elapsed_refL)
{
    double numerator = 4.0/3.0*k*k*k - k*k*(2.0*m+2.0*n+3.0) + k*m*(4.0*n+3.0);
    return (numerator)/(elapsed_refL*1.0e+9);
}
double computeKPerform_refL(double m, double n, double elapsed_refL)
{
    double kr = 2.0*m*n*n + n*n - n;
    kr /= 2.0;
    double ut = 2.0*m*m*m + 3.0*m*m + 6.0*m*n*n + 6.0*m*n - 5.0*m + 2.0*n*n*n - 12.0*n*n + 10.0*n;
    ut /= 6;
    double numerator = kr + ut;
    return (numerator)/(elapsed_refL*1.0e+9);
}

int main(int argc, char *argv[]) {
    // integer variables
    int info, lda, ldq, m, n, k, lwork, nb, i, j;
    int workQuery = -1;
    // double variables
    double *A, *Q, *As, *tau, *work, *T=NULL;
    double normA, tmp;
    double perform_refL, elapsed_refL, norm_orth_1, norm_repres_1;
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
    n = 20;
    k = -1;
    nb = 3; // Choose a default nb value that is NOT a factor of k
    while (k % nb == 0) {
        nb++;
    }
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
        if( strcmp( *(argv + i), "-nb") == 0) {
            nb  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-t") == 0) {
            timesOnly = true;
        }
    }

    if( lda < 0 ) lda = m;
    if( ldq < 0 ) ldq = m;
    // While our functionality works for when k < n,we are constructing a different
    // matrix Q, so to keep the comparisons reasonable, we will enforce k=n
    k = n;

    // allocate memory for the matrices and vectors that we need to use
    A =   (double *) malloc(lda * k * sizeof(double));
    As =  (double *) malloc(lda * k * sizeof(double));
    Q =   (double *) malloc(ldq * k * sizeof(double));
    tau = (double *) malloc(k * sizeof(double));
    T =   (double *) malloc(n * n * sizeof(double));

    // Print to the user what we are doing along with any arguments that are used
    char *source = SOURCE;
    if(!timesOnly) {
        printf("dgeqrf dorg2r %s | m = %4d, n = %d, k = %4d, lda = %4d, ldq = %4d\n", source, m, n, k, lda, ldq);
    }

    // Generate A as a random matrix
    for (i = 0; i < lda * k; i++)
        A[i] = (double) rand() / (double) (RAND_MAX) - 0.5e+00;
    // Store random data inside Q to ensure that we do not assume anything about Q
    for (i = 0; i < ldq * k; i++)
        Q[i] = (double) rand() / (double) (RAND_MAX) - 0.5e+00;

    // Store a copy of A inside As
    dlacpy_(&aChar, &m, &k, A, &lda, As, &lda, dummy);
    // Find the norm of A for use in later accuracy computations
    normA = 0.0;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            tmp = A[i + j*m];
            normA += tmp*tmp;
        }
    }
    normA = sqrt(normA);
    // Create the work array to do workspace queries
    work = (double *) malloc(sizeof(double));
    // Determine how much workspace is needed for our operations
    // Check dgeqrf first
    dgeqrf_(&m, &k, A, &lda, tau, work, &workQuery, &info );
    lwork = work[0];

    // reallocate work to be of the right size
    work = (double *) realloc(work, lwork * sizeof(double));

    // Call dgeqrf first
    dgeqrf_(&m, &k, A, &lda, tau, work, &lwork, &info);

    // Copy A into Q for use with dorgkr
    dlacpy_(&aChar, &m, &k, A, &lda, Q, &ldq, dummy);

    // Take the current time for use with timing dorg2r
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    // Construct the T matrix associated with the householder reflectors
    my_dlarft_ut_(&m, &k, Q, &ldq, tau, Q, &ldq);
    //dlarft_(&fChar, &cChar, &m, &k, Q, &ldq, tau, Q, &ldq, dummy, dummy);

    // Copy the upper triangular part of T into the upper triangular part of Q (Overwriting R)
    //dlacpy_(&uChar, &n, &n, T, &n, Q, &ldq, dummy);

    // Call my_dorgkr
    my_dorgkr_(&m, &k, Q, &ldq);
    //dorg2r_(&m, &n, &k, Q, &ldq, tau, work, &info);

    // Determine how much time has taken
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    // Compute the error information
    // reallocate work to be of the right size for testing orthogonality of Q.
    work = (double *)realloc(work, k * k * sizeof(double));

    // Set work to be I
    dlaset_(&aChar, &k, &k, &zero, &one, work, &k, dummy);
    // Compute work = Q**T * Q - I
    dsyrk_(&uChar, &tChar, &k, &m, &one, Q, &ldq, &negOne, work, &k);
    // Compute the norm of work
    for (i=0; i < k; i++) {
        for (j=0; j < k; j++) {
            tmp = work[i + j*k] * work[i + j*k];
            norm_orth_1 += tmp;
        }
    }
    norm_orth_1 = sqrt(norm_orth_1);

    // reallocate work to be able to hold Q
    work = (double *)realloc(work, m * k * sizeof(double));
    // Copy Q into work
    dlacpy_(&aChar, &m, &k, Q, &ldq, work, &m, dummy);

    // Compute work = work * R
    dtrmm_(&rChar, &uChar, &nChar, &nChar, &m, &k, &one, A, &lda, work, &m, dummy, dummy, dummy, dummy);
    
    // Compute work = work - A
    for (i = 0; i < m; i++) {
        for (j = 0; j < k; j++) {
            tmp = work[i +j*m] - As[i+j*lda];
            norm_repres_1 += tmp * tmp;
        }
    }
    // Compute ||A - QR||_F
    norm_repres_1 = sqrt(norm_repres_1);
    // Compute ||A - Q*R||_F / ||A||_F
    norm_repres_1 /= normA;

    // Compute 'performance'
    // This is the performance of our algorithm. It is a much different algorithm than 
    perform_refL = computeKPerform_refL((double) m, (double) n, elapsed_refL);
    if(!timesOnly) {
        printf("my_dorgkr: perf = %f time = %f repres = %5.1e ortho = %5.1e", perform_refL, elapsed_refL, norm_repres_1, norm_orth_1);
    } else {
        printf("k: %10.10e:%10.10e", elapsed_refL, perform_refL);
    }
    printf("\n");

    // Now do the regular dorg2r to compare time.
    // Since we are not testing the accuracy of this method, we only need to call dorg2r correctly
    work = (double *) realloc(work, 1 * sizeof(double));

    lwork = -1;
    //dorgqr_(&m, &n, &k, As, &lda, tau, work, &lwork, &info);
    lwork = n;
    work = (double *) realloc(work, lwork * sizeof(double));
    // Take the current time for use with timing dorg2r
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    // Call dorg2r
    //dorgqr_(&m, &n, &k, As, &lda, tau, work, &lwork, &info);
    dorg2r_(&m, &k, &k, As, &lda, tau, work, &info);

    // Determine how much time has taken
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    perform_refL = computePerform_refL((double) m, (double) n, (double) k, elapsed_refL);
    if(!timesOnly) {
        printf("aocl dorg2r: perf = %f time = %f", perform_refL, elapsed_refL);
    } else {
        printf("a: %10.10e:%10.10e", elapsed_refL, perform_refL);
    }
    printf("\n");
    //printf("dorg2r: time = %f repres = %5.1e ortho = %5.1e\n", elapsed_refL, norm_repres_1, norm_orth_1);
    // Free memory before terminating 
    // Take the current time for use with timing dorg2r
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    // Call dorg2r
    //dorgqr_(&m, &n, &k, As, &lda, tau, work, &lwork, &info);
    ref_dorg2r_(&m, &k, &k, As, &lda, tau, work, &info);

    // Determine how much time has taken
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    perform_refL = computePerform_refL((double) m, (double) n, (double) k, elapsed_refL);
    if(!timesOnly) {
        printf("reference dorg2r: perf = %f time = %f", perform_refL, elapsed_refL);
    } else {
        printf("r: %10.10e:%10.10e", elapsed_refL, perform_refL);
    }
    printf("\n");
    free(Q);
    free(As);
    free(A);
    free(tau);
    free(work);
    free(T);
}
