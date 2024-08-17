// This file is responsible for timing my_dlarfb vs dlarfb
//
// TODO:
// dgeqrf -> dlarft -> {my_,}dlarfb -> A'
// compare outputs of dlarfb (compare different A')
// apply dorg2r and then compare results
//
//
// timing: just calls to {my_,}dlarfb
// TODO:
// accuracy: both after {my_,}dlarfb and dorg2r
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

void checkMemory(void *ptr, char *name, size_t size)
{
    if (ptr == NULL) {
        fprintf(stderr, "Insufficient Memory. Attempted to allocate %ld bytes for array %s\n", size, name);
        exit(1);
    }
}

int main(int argc, char **argv)
{
    // Local params
    int info, m, n, k, lwork;
    double *Q, *Qs, *C, *Cs, *T, *workMat, *tau, *work=NULL;
    double normA;
    double elapsed_refL, perform_refL;
    struct timeval tp;
    int i, j;

    char aChar = 'A';
    char cChar = 'C';
    char fChar = 'F';
    char lChar = 'L';
    char nChar = 'N';
    char uChar = 'U';
    int dummy = 0;

    int negOne = -1;
    double one = 1.0;
    double zero = 0.0;
    double dNegOne = -1.0;

    // Init random seed
    //srand(0);

    // Default parameters to ensure that we hit 
    // the blocked code
    m = 30;
    n = 20;
    k = n/2 + 1;
    // Flag that helps facilitate testing with the driver.c file
    bool errorsOnly = false;
    // Flag that helps facilitate timing with the time.sh file
    bool timesOnly = false;

    for(i = 1; i < argc; ++i){
        if( strcmp( *(argv + i), "-m") == 0) {
            m  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-n") == 0) {
            n  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-k") == 0) {
            k  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-t") == 0) {
            timesOnly  = true;
        }
    }

    // Guard against making k too big on accident
    if( k > n ) k = (n == 1) ? n : n/2 + 1;

    if (!errorsOnly && !timesOnly) {
        printf("dgeqrf dlarft dlarfb LAPACK | ");
        printf("m = %4d, ",    m);
        printf("n = %4d, ",    n);
        printf("k = %4d\n",    k);
        //printf("             ");
    }

    Q  = (double *) calloc(m * n, sizeof(double));
    Qs  = (double *) calloc(m * n, sizeof(double));
    T  = (double *) calloc(k * k, sizeof(double));
    C  = (double *) calloc(m * n, sizeof(double));
    Cs  = (double *) calloc(m * n, sizeof(double));
    checkMemory(Q, "Q", m * n * sizeof(double));
    checkMemory(Qs, "Qs", m * n *  sizeof(double));
    checkMemory(T, "T", k * k * sizeof(double));
    checkMemory(C, "C", m * n * sizeof(double));
    checkMemory(Cs, "Cs", m * n * sizeof(double));
    // Part of the call to dlarfb requires a workspace
    // may be able to refactor it out, but not really worth the effort on a timing test file
    workMat = (double *) calloc(m * k, sizeof(double));
    checkMemory(workMat, "workMat", m * k * sizeof(double));

    for(i = 0; i < m * n; ++i)
        *(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

    // Populate C with a random matrix for use with computing 
    // time with dlarfb.
    for(i = 0; i < m * n; ++i)
        *(C + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
    dlacpy_(&aChar, &m, &n, C, &m, Cs, &m, dummy);

    tau = (double *) malloc( n * sizeof(double));
    checkMemory(tau, "tau", n * sizeof(double));

    work = (double *) malloc( 1 * sizeof(double));
    checkMemory(work, "work", sizeof(double));
    dgeqrf_(&m, &n, Q, &m, tau, work, &negOne, &info);
    //LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, k, A, lda, tau, work, -1 ); 
    lwork = ((int) work[0]);
    free( work );
    work = (double *) malloc( lwork * sizeof(double));
    checkMemory(work, "work", lwork * sizeof(double));
    
    // Compute the householder reflectors for A
    dgeqrf_(&m, &n, Q, &m, tau, work, &lwork, &info);
    // Store these reflectors
    dlacpy_(&aChar, &m, &n, Q, &m, Qs, &m, dummy);
    // Construct the associated T factor of these reflectors
    dlarft_(&fChar, &cChar, &m, &k, Q, &m, tau, T, &k, dummy, dummy);
    // Now, we create the timing struct
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    // Call existing dlarfb
    dlarfb_(&lChar, &nChar, &fChar, &cChar, &m, &n, &k, Q, &m, T, &k, C, &m, workMat, &m);
    
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    // print out time of optimized dlarfb
    printf("OPT: %10.10e\n", elapsed_refL);
    // copy everything back to how it was for Q and C (T should not be touched and we are trusting
    // the optimized implementers to do this properly)
    dlacpy_(&aChar, &m, &n, Qs, &m, Q, &m, dummy);
    dlacpy_(&aChar, &m, &n, Cs, &m, C, &m, dummy);
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    dlarfb_(&lChar, &nChar, &fChar, &cChar, &m, &n, &k, Q, &m, T, &k, C, &m, workMat, &m);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    printf("REF:  %10.10e\n", elapsed_refL);
    // copy everything back to how it was for Q and C (T should not be touched and we are trusting
    // the optimized implementers to do this properly)
    dlacpy_(&aChar, &m, &n, Qs, &m, Q, &m, dummy);
    dlacpy_(&aChar, &m, &n, Cs, &m, C, &m, dummy);
    dlacpy_(&uChar, &k, &k, T, &k, Q, &m, dummy);
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    my_dlarfb_(&m, &n, &k, Q, &m, C, &m);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    printf("MY:  %10.10e\n", elapsed_refL);
}
