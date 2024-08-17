#include <stdio.h>
#include<stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

double myNorm(int m, int n, double *A, int lda)
{
    double ret = 0.0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double tmp = A[i + j * lda];
            ret += tmp*tmp;
        }
    }
    return sqrt(ret);
}

int main(int argc, char **argv) {

    // Local params
    int info, lda, ldq, m, n, k, lwork, nb;
    double elapsed_refL, perform_refL;
    struct timeval tp;
    int i, j, version;

    char aChar = 'A';
    char fChar = 'F';
    char uChar = 'U';
    char tChar = 'T';
    char rChar = 'R';
    char nChar = 'N';
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
    nb = 3; // Choose a default nb value that is NOT a factor of k
    lda = -1;
    ldq = -1;
    version = -1; // Default to the system version (vendor usually)
    // Flag that helps facilitate testing with the driver.c file
    bool errorsOnly = false;
    // Flag that helps facilitate timing with the time.sh file
    bool timesOnly = false;

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
        if( strcmp( *(argv + i), "-k") == 0) {
            k  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-e") == 0) {
            errorsOnly  = true;
        }
        if( strcmp( *(argv + i), "-t") == 0) {
            timesOnly  = true;
        }
        if( strcmp( *(argv + i), "-v") == 0) {
            version = atoi( *(argv + i + 1) );
            i++;
        }
    }

    if( lda < 0 ) lda = m;
    if( ldq < 0 ) ldq = m;
    if( k > n) k = n;

    if (!errorsOnly && !timesOnly) {
        printf("dgeqrf dorgqr | ");
        printf("version = %4d, ", version);
        printf("m = %4d, ",    m);
        printf("n = %4d, ",    n);
        printf("k = %4d, ",    k);
        printf("nb = %4d, ",    nb);
        printf("lda = %4d, ",lda);
        printf("ldq = %4d, ",ldq);
        printf("\n");
    }

    double *A  = (double *) calloc(lda * k,sizeof(double));
    double *As = (double *) calloc(lda * k, sizeof(double));
    double *Q  = (double *) calloc(ldq * n, sizeof(double));

    for(i = 0; i < lda * k; ++i)
        *(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

    for(i = 0; i < ldq * n; ++i)
        *(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

    dlacpy_(&aChar, &m, &k, A, &lda, As, &lda, dummy);
    double normA = myNorm(m, k, A, lda);
    //info  = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, k, A, lda, As, lda );
    //normA = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, k, A, lda, work );

    double *tau = (double *) malloc( k * sizeof(double));

    double *work = (double *) malloc( 1 * sizeof(double));
    dgeqrf_(&m, &k, A, &lda, tau, work, &negOne, &info);
    //LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, k, A, lda, tau, work, -1 ); 
    lwork = ((int) work[0]);
    // determining what function to call based on the value of 'version'
    switch(version) 
    {
        case 0:
            my_dorgqr_v0_(&m, &n, &k, Q, &ldq, tau, work, &negOne, &info);
            break;
        case 1:
            my_dorgqr_v1_(&m, &n, &k, Q, &ldq, tau, work, &negOne, &info);
            break;
        case 2:
            my_dorgqr_v2_(&m, &n, &k, Q, &ldq, tau, work, &negOne, &info);
            break;
        case 3:
            my_dorgqr_v3_(&m, &n, &k, Q, &ldq, tau, work, &negOne, &info);
            break;
        case 4:
            my_dorgqr_v4_(&m, &n, &k, Q, &ldq, tau, work, &negOne, &info);
            break;
        case 5:
            my_dorgqr_v5_(&m, &n, &k, Q, &ldq, tau, work, &negOne, &info);
            break;
        case 6:
            my_dorgqr_v6_(&m, &n, &k, Q, &ldq, tau, work, &negOne, &info);
            break;
        case 7:
            my_dorgqr_v7_(&m, &n, &k, Q, &ldq, tau, work, &negOne, &info);
            break;
        case 8:
            my_dorgqr_v8_(&m, &n, &k, Q, &ldq, tau, work, &negOne, &info);
            break;
        case 9:
            my_dorgqr_v9_(&m, &n, &k, Q, &ldq, tau, work, &negOne, &info);
            break;
        case 10:
            my_dorgqr_v10_(&m, &n, &k, Q, &ldq, tau, work, &negOne, &info);
            break;
        case 11:
            my_dorgqr_v11_(&m, &n, &k, Q, &ldq, tau, work, &negOne, &info);
            break;
        default:
            dorgqr_(&m, &n, &k, Q,&ldq, tau, work, &negOne, &info);
    }
    //LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, k, A, lda, tau, work, -1 );
    if (lwork < ((int) work[0])) lwork = ((int) work[0]); 
    free( work );
    if (lwork < n) lwork = n;
    work = (double *) malloc( lwork * sizeof(double));

    dgeqrf_(&m, &k, A, &lda, tau, work, &lwork, &info);
    //LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, k, A, lda, tau, work, lwork ); 
    dlacpy_(&aChar, &m, &k, A, &lda, Q, &ldq, dummy);
    //LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, k, A, lda, Q, ldq );
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    //LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, k, Q, ldq, tau, work, lwork );
    // Directly calling the fortran function dorgqr
    // Directly calling my fortran function my_dorgqr
    // determining what function to call based on the value of 'version'
    switch(version) 
    {
        case 0:
            my_dorgqr_v0_(&m, &n, &k, Q, &ldq, tau, work, &lwork, &info);
            break;
        case 1:
            my_dorgqr_v1_(&m, &n, &k, Q, &ldq, tau, work, &lwork, &info);
            break;
        case 2:
            my_dorgqr_v2_(&m, &n, &k, Q, &ldq, tau, work, &lwork, &info);
            break;
        case 3:
            my_dorgqr_v3_(&m, &n, &k, Q, &ldq, tau, work, &lwork, &info);
            break;
        case 4:
            my_dorgqr_v4_(&m, &n, &k, Q, &ldq, tau, work, &lwork, &info);
            break;
        case 5:
            my_dorgqr_v5_(&m, &n, &k, Q, &ldq, tau, work, &lwork, &info);
            break;
        case 6:
            my_dorgqr_v6_(&m, &n, &k, Q, &ldq, tau, work, &lwork, &info);
            break;
        case 7:
            my_dorgqr_v7_(&m, &n, &k, Q, &ldq, tau, work, &lwork, &info);
            break;
        case 8:
            my_dorgqr_v8_(&m, &n, &k, Q, &ldq, tau, work, &lwork, &info);
            break;
        case 9:
            my_dorgqr_v9_(&m, &n, &k, Q, &ldq, tau, work, &lwork, &info);
            break;
        case 10:
            my_dlarft_ut_(&m, &k, Q, &ldq, tau, Q, &ldq); 
            my_dorgqr_v10_(&m, &n, &k, Q, &ldq, tau, work, &lwork, &info);
            break;
        case 11:
            my_dlarft_ut_(&m, &k, Q, &ldq, tau, Q, &ldq); 
            my_dorgqr_v11_(&m, &n, &k, Q, &ldq, tau, work, &lwork, &info);
            break;
        default:
            dorgqr_(&m, &n, &k, Q,&ldq, tau, work, &lwork, &info);
    }

    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    free( work );
    // fix later. This was from timing everything
    perform_refL = (4.0e+00*(double) m*(double)n*(double)k -2.0e+00*(double)m*(double)k*(double)k - 2.0e+00*(double)n*(double)k*(double)k + 4.0e+00/3.0e+00 *(double)k*(double)k*(double)k)/elapsed_refL /1.0e+9;
    
    double norm_orth_1, norm_repres_1;

    work  = (double *) malloc(n * n * sizeof(double));
    dlaset_(&aChar, &n, &n, &zero, &one, work, &n, dummy);
    //info  = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (1e+00), work, n );
    //bl1_dsyrk_blas(CblasColMajor, CblasUpper, CblasTrans, n, m, 1.0e+00, Q, ldq, -1.0e+00, work, n);
    dsyrk_(&uChar, &tChar, &n, &m, &one, Q, &ldq, &dNegOne, work, &n, dummy, dummy);
    //cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, m, 1.0e+00, Q, ldq, -1.0e+00, work, n );
    norm_orth_1 = myNorm(n, n, work, n);
    //norm_orth_1 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, n, work, n, NULL );
    free( work );

    work  = (double *) malloc(m * k * sizeof(double));
    dlacpy_(&aChar, &m, &k, Q, &ldq, work, &m, dummy);
    //info  = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, k, Q, ldq, work, m );
    dtrmm_(&rChar, &uChar, &nChar, &nChar, &m, &k, &one, A, &lda, work, &m, dummy, dummy, dummy, dummy);
    //cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, k, (1.0e+00), A, lda, work, m );
    for(i = 0; i < m; ++i)
        for(j = 0; j < k; ++j)
            work[ i+j*m ] -= As[ i+j*lda ];
    norm_repres_1 = myNorm(m, k, work, m);
    //norm_repres_1 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, k, work, m, NULL );
    norm_repres_1 = norm_repres_1 / normA;
    free( work );

    if (!errorsOnly && !timesOnly) {
        printf("| time = %f   GFlop/sec = %f", elapsed_refL, perform_refL);
        printf("| repres  = %5.1e    ortho = %5.1e ", norm_repres_1, norm_orth_1);
    } else if (errorsOnly && !timesOnly) {
        printf("%10.10e %10.10e", norm_repres_1, norm_orth_1);
    } else {
        printf("%10.10e:%10.10e", elapsed_refL, perform_refL);
    }

    printf("\n");

    free( Q );
    free( A );
    free( As );
    free( tau );

    return 0;
}
