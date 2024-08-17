#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <cblas.h>
#include <lapacke.h>

//Since this isn't defined in regular C
int min(int a, int b) {
    return ((a < b) ? a : b);
}

int main(int argc, char *argv[]) {
    int info, lda, ldq, m, n, k, lwork, nb;
    double *A, *U, *As, *s, *work=NULL;
    double normA;
    double elapsed_refL, perform_refL;
    struct timeval tp;
    int i, j;
    // In order to time the SVD we will first construct a matrix then compute its SVD while only
    // asking for the singular values and left singular vectors

    // Default parameters to ensure that we hit 
    // the blocked code
    m = 30;
    n = 20;
    lda = -1;
    for(i = 1; i < argc; ++i){
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
    }

    if( lda < 0 ) lda = m;

    A = (double *) calloc(lda * n, sizeof(double));
    U = (double *) calloc(lda * m, sizeof(double));
    As = (double *) calloc(lda * n, sizeof(double));
    s = (double *) calloc(min(m,n), sizeof(double));
    for(i = 0; i < lda * n; ++i)
        *(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
    // Save a copy of A in As
    info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, As, lda);
    //----------------------------------------------
    // unmodified DORGQR
    //----------------------------------------------
    char jobU = 'A'; // Grab the left singular vectors
    char jobVt= 'N'; // Do not compute the right singular vectors
    // Start the timer
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    // Workspace computation.
    lwork = -1;
    work = (double *) malloc(2 * sizeof(double));
    size_t dummy;
    dgesvd_(&jobU, &jobVt, &m, &n, A, &lda, s, U, &lda, A, &lda, work, &lwork, &info, dummy, dummy);
    lwork = work[0];
    free(work);
    // Actually compute the SVD
    work = (double *) malloc(lwork * sizeof(double));
    dgesvd_(&jobU, &jobVt, &m, &n, A, &lda, s, U, &lda, A, &lda, work, &lwork, &info, dummy, dummy);
    // Stop the timer
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    printf("m = %d, n = %d\n", m, n);
    printf("Time for unmodified DORGQR:\t%lf\n", elapsed_refL);
    //----------------------------------------------
    // MY_DORGQR
    //----------------------------------------------
    free(work);
    // Start the timer
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    // Workspace computation.
    lwork = -1;
    work = (double *) malloc(2*sizeof(double));
    my_dgesvd_(&jobU, &jobVt, &m, &n, As, &lda, s, U, &lda, A, &lda, work, &lwork, &info);
    lwork = work[0];
    free(work);
    // Actually compute the SVD
    work = (double *) malloc(lwork * sizeof(double));
    my_dgesvd_(&jobU, &jobVt, &m, &n, As, &lda, s, U, &lda, A, &lda, work, &lwork, &info);
    // Stop the timer
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    printf("Time for modified DORGQR:\t%lf\n", elapsed_refL);
    //----------------------------------------------
    // Potential TODO: ensure the singular vectors computed 
    //  using each method is the same. 
    // They should be the same because of results from
    //  
    //----------------------------------------------
    // Free memory
    free(A);
    free(U);
    free(As);
    free(s);

}
