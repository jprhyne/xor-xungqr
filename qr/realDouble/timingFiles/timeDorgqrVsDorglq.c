#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "../test/perf_calc.h"
int main(int argc, char **argv) {
    // Local params
    int info, m, n, k, lwork;
    size_t mS, nS, kS, iS;
    double *A, *Q, *Qt, *tau1, *tau2, *workMat;
    double normA, tmpVal;
    double timeQR, perfQR, timeLQ, perfLQ;
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

    // Default parameters to ensure that we hit 
    // the blocked code
    m = 30;
    n = 20;
    k = n/2 + 1;
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
    }
    // Enforce m \geq n \geq k
    if (m < n) {
        printf("m must be equal to or larger than n\n");
        return 1;
    }
    if (n < k) {
        printf("n must be equal to or larger than k\n");
        return 2;
    }
    mS = (size_t) m;
    nS = (size_t) n;
    kS = (size_t) k;

    // Guard against making k too big on accident
    if( k > n ) k = (n == 1) ? n : n/2 + 1;

    // Allocate our matrices
    A = (double *) calloc(mS*kS, sizeof(double));
    Q = (double *) calloc(mS*nS, sizeof(double));
    Qt= (double *) calloc(nS*mS, sizeof(double));
    tau1 = (double *) calloc(nS, sizeof(double));
    tau2 = (double *) calloc(nS, sizeof(double));

    // Generate the matrix A to be m by k (where m \geq n\geq k)
    for(iS = 0; iS < mS * kS; ++iS)
        *(A + iS) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

    // Copy A into Q and A**T into Qt
    for(i = 0; i < m; ++i) {
        for (j=0; j < k; ++j) {
            Q[i + j*m] = A[i + j*m];
            Qt[j + i*n] = A[i + j*m];
        }
    }
    
    // Determine the workspace size requirements for the entire path
    double workQuery[1] = {0};
    dgeqrf_(&m, &n, Q, &m, tau1, workQuery, &negOne, &info);
    lwork = ((int) workQuery[0]);
    dgelqf_(&n, &m, Qt, &n, tau2, workQuery, &negOne, &info);
    lwork = (lwork < ((int) workQuery[0])) ? ((int) workQuery[0]) : lwork;
    dorgqr_lap_(&m, &n, &k, Q, &m, tau1, workQuery, &negOne, &info);
    lwork = (lwork < ((int) workQuery[0])) ? ((int) workQuery[0]) : lwork;
    dorglq_lap_(&n, &m, &k, Qt, &n, tau2, workQuery, &negOne, &info);
    lwork = (lwork < ((int) workQuery[0])) ? ((int) workQuery[0]) : lwork;

    // Allocate our workspace to be the needed size
    workMat = (double *) malloc(lwork * sizeof(double));
    // factorize A and A**T
    dgeqrf_(&m, &n, Q,  &m, tau1, workMat, &lwork, &info);
    dgelqf_(&n, &m, Qt, &n, tau2, workMat, &lwork, &info);
    // Compute A = QR
    gettimeofday(&tp, NULL);
    timeQR=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    dorgqr_lap_(&m, &n, &k, Q, &m, tau1, workMat, &lwork, &info);
    gettimeofday(&tp, NULL);
    timeQR+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    // Compute A**T = LQ
    gettimeofday(&tp, NULL);
    timeLQ=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    dorglq_lap_(&n, &m, &k, Qt, &n, tau1, workMat, &lwork, &info);
    gettimeofday(&tp, NULL);
    timeLQ+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    // Compute our performance values
    perfQR = computeDorgqxPerf((double) m, (double) n, (double) k, timeQR);
    perfLQ = computeDorgxqPerf((double) n, (double) m, (double) k, timeLQ);

    // Print to console
    printf("qr:New:%10.10e\n",perfQR);
    printf("lq:New:%10.10e\n",perfLQ);

    // Copy A into Q and A**T into Qt
    for(i = 0; i < m; ++i) {
        for (j=0; j < k; ++j) {
            Q[i + j*m] = A[i + j*m];
            Qt[j + i*n] = A[i + j*m];
        }
    }
    dgeqrf_(&m, &n, Q, &m, tau1, workQuery, &negOne, &info);
    lwork = ((int) workQuery[0]);
    dgelqf_(&n, &m, Qt, &n, tau2, workQuery, &negOne, &info);
    lwork = (lwork < ((int) workQuery[0])) ? ((int) workQuery[0]) : lwork;
    dorgqr_ref_(&m, &n, &k, Q, &m, tau1, workQuery, &negOne, &info);
    lwork = (lwork < ((int) workQuery[0])) ? ((int) workQuery[0]) : lwork;
    dorglq_ref_(&n, &m, &k, Qt, &n, tau2, workQuery, &negOne, &info);
    lwork = (lwork < ((int) workQuery[0])) ? ((int) workQuery[0]) : lwork;
    free(workMat);

    // Allocate our workspace to be the needed size
    workMat = (double *) malloc(lwork * sizeof(double));
    // factorize A and A**T
    dgeqrf_(&m, &n, Q,  &m, tau1, workMat, &lwork, &info);
    dgelqf_(&n, &m, Qt, &n, tau2, workMat, &lwork, &info);
    // Compute A = QR
    gettimeofday(&tp, NULL);
    timeQR=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    dorgqr_ref_(&m, &n, &k, Q, &m, tau1, workMat, &lwork, &info);
    gettimeofday(&tp, NULL);
    timeQR+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    // Compute A**T = LQ
    gettimeofday(&tp, NULL);
    timeLQ=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    dorglq_ref_(&n, &m, &k, Qt, &n, tau1, workMat, &lwork, &info);
    gettimeofday(&tp, NULL);
    timeLQ+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    // Compute our performance values
    perfQR = computeDorgqxPerf((double) m, (double) n, (double) k, timeQR);
    perfLQ = computeDorgxqPerf((double) n, (double) m, (double) k, timeLQ);

    // Print to console
    printf("qr:Old:%10.10e\n",perfQR);
    printf("lq:Old:%10.10e\n",perfLQ);

    free(A);
    free(Q);
    free(Qt);
    free(tau1);
    free(tau2);
    free(workMat);
}
