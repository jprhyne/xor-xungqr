#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <complex.h>
#include <sys/time.h>
void setAf(int m, int n, float *A) {
    for (size_t i = 0; i < m*n; ++i) {
        A[i] = (float)rand() / (float)(RAND_MAX) - 0.5e+00;
    }
}
void setAd(int m, int n, double *A) {
    for (size_t i = 0; i < m*n; ++i) {
        A[i] = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
    }
}
void setAfc(int m, int n, float complex *A) {
    for (size_t i = 0; i < m*n; ++i) {
        float real = (float)rand() / (float)(RAND_MAX) - 0.5e+00;
        float imag = (float)rand() / (float)(RAND_MAX) - 0.5e+00;
        A[i] = real + imag * I;
    }
}
void setAdc(int m, int n, double complex *A) {
    for (size_t i = 0; i < m*n; ++i) {
        double real = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
        double imag = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
        A[i] = real + imag * I;
    }
}
/*
 * On entry, timeVals must be a 2d double array of dimension timeVals[4][3]
 * timeVals has the following entries on exit:
 * timeVals[0][0]: qr AOCL time
 * timeVals[0][1]: qr reference time
 * timeVals[0][2]: qr my time
 * timeVals[1][0]: ql AOCL time
 * timeVals[1][1]: ql reference time
 * timeVals[1][2]: ql my time
 * timeVals[2][0]: rq AOCL time
 * timeVals[2][1]: rq reference time
 * timeVals[2][2]: rq my time
 * timeVals[3][0]: lq AOCL time
 * timeVals[3][1]: lq reference time
 * timeVals[3][2]: lq my time
 */
void timeReal(int m, int n, int k, double timeVals[4][3]) {
    struct timeval tp;

    double elapsed_refL;

    float *A = (float *) malloc(sizeof(float)*m*k);
    float *Q = (float *) malloc(sizeof(float)*m*n);
    float *tau = (float *) malloc(sizeof(float) * (m < n ? n : n)); 
    float *work = (float *) malloc(sizeof(float));

    int neg_one = -1;
    int lwork = -1;
    int info = 0;
    // set all of timeVals to -1 (to ensure everything is set on exit)
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            timeVals[i][j] = -1;
        }
    }
    // workspace queries for all factorization routines
    sgeqrf_(&m, &k, A, &m, tau, work, &neg_one, &info);
    lwork = (int) work[0];
    sgeqlf_(&m, &k, A, &m, tau, work, &neg_one, &info);
    if (lwork < (int) work[0]) lwork = (int) work[0];
    sgerqf_(&k, &m, A, &k, tau, work, &neg_one, &info);
    if (lwork < (int) work[0]) lwork = (int) work[0];
    sgelqf_(&k, &m, A, &k, tau, work, &neg_one, &info);
    if (lwork < (int) work[0]) lwork = (int) work[0];
    // free the workspace
    free(work);
    // allocate the workspace to be as large as we need
    work = (float *) malloc(sizeof(float) * lwork);
    //----------------------------------------------------------------------------------------------
    // qr
    //----------------------------------------------------------------------------------------------
    setAf(m, k, A);
    sgeqrf_(&m, &k, A, &m, tau, work, &lwork, &info);
    //----------------------------------------------------------------------------------------------
    // AOCL
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    sorgqr_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[0][0] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // Old reference
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    sorgqr_ref_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[0][1] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // My version
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    sorgqr_new_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[0][2] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // ql
    //----------------------------------------------------------------------------------------------
    setAf(m, k, A);
    sgeqlf_(&m, &k, A, &m, tau, work, &lwork, &info);
    //----------------------------------------------------------------------------------------------
    // AOCL
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    sorgql_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[1][0] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // Old reference
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    sorgql_ref_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[1][1] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // My version
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    sorgql_new_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[1][2] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // rq
    //----------------------------------------------------------------------------------------------
    setAf(m, k, A);
    sgerqf_(&k, &m, A, &k, tau, work, &lwork, &info);
    //----------------------------------------------------------------------------------------------
    // AOCL
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    sorgrq_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[2][0] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // Old reference
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    sorgrq_ref_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[2][1] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // My version
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    sorgrq_new_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[2][2] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // lq
    //----------------------------------------------------------------------------------------------
    setAf(m, k, A);
    sgelqf_(&k, &m, A, &k, tau, work, &lwork, &info);
    //----------------------------------------------------------------------------------------------
    // AOCL
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    sorglq_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[3][0] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // Old reference
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    sorglq_ref_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[3][1] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // My version
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    sorglq_new_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[3][2] = elapsed_refL;
}

void timeDouble(int m, int n, int k, double timeVals[4][3]) {
    struct timeval tp;

    double elapsed_refL;

    double *A = (double *) malloc(sizeof(double)*m*k);
    double *Q = (double *) malloc(sizeof(double)*m*n);
    double *tau = (double *) malloc(sizeof(double) * (m < n ? n : n)); 
    double *work = (double *) malloc(sizeof(double));

    int neg_one = -1;
    int lwork = -1;
    int info = 0;
    // set all of timeVals to -1 (to ensure everything is set on exit)
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            timeVals[i][j] = -1;
        }
    }
    // workspace queries for all factorization routines
    dgeqrf_(&m, &k, A, &m, tau, work, &neg_one, &info);
    lwork = (int) work[0];
    dgeqlf_(&m, &k, A, &m, tau, work, &neg_one, &info);
    if (lwork < (int) work[0]) lwork = (int) work[0];
    dgerqf_(&k, &m, A, &k, tau, work, &neg_one, &info);
    if (lwork < (int) work[0]) lwork = (int) work[0];
    dgelqf_(&k, &m, A, &k, tau, work, &neg_one, &info);
    if (lwork < (int) work[0]) lwork = (int) work[0];
    // free the workspace
    free(work);
    // allocate the workspace to be as large as we need
    work = (double *) malloc(sizeof(double) * lwork);
    //----------------------------------------------------------------------------------------------
    // qr
    //----------------------------------------------------------------------------------------------
    setAd(m, k, A);
    dgeqrf_(&m, &k, A, &m, tau, work, &lwork, &info);
    //----------------------------------------------------------------------------------------------
    // AOCL
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    dorgqr_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[0][0] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // Old reference
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    dorgqr_ref_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[0][1] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // My version
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    dorgqr_new_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[0][2] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // ql
    //----------------------------------------------------------------------------------------------
    setAd(m, k, A);
    dgeqlf_(&m, &k, A, &m, tau, work, &lwork, &info);
    //----------------------------------------------------------------------------------------------
    // AOCL
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    dorgql_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[1][0] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // Old reference
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    dorgql_ref_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[1][1] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // My version
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    dorgql_new_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[1][2] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // rq
    //----------------------------------------------------------------------------------------------
    setAd(m, k, A);
    dgerqf_(&k, &m, A, &k, tau, work, &lwork, &info);
    //----------------------------------------------------------------------------------------------
    // AOCL
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    dorgrq_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[2][0] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // Old reference
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    dorgrq_ref_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[2][1] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // My version
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    dorgrq_new_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[2][2] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // lq
    //----------------------------------------------------------------------------------------------
    setAd(m, k, A);
    dgelqf_(&n, &m, A, &n, tau, work, &lwork, &info);
    //----------------------------------------------------------------------------------------------
    // AOCL
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    dorglq_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[3][0] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // Old reference
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    dorglq_ref_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[3][1] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // My version
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    dorglq_new_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[3][2] = elapsed_refL;
}

void timeCReal(int m, int n, int k, double timeVals[4][3]) {
    struct timeval tp;

    double elapsed_refL;

    float complex *A = (float complex *) malloc(sizeof(float complex)*m*k);
    float complex *Q = (float complex *) malloc(sizeof(float complex)*m*n);
    float complex *tau = (float complex *) malloc(sizeof(float complex) * (m < n ? n : n)); 
    float complex *work = (float complex *) malloc(sizeof(float complex));

    int neg_one = -1;
    int lwork = -1;
    int info = 0;
    // set all of timeVals to -1 (to ensure everything is set on exit)
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            timeVals[i][j] = -1;
        }
    }
    // workspace queries for all factorization routines
    cgeqrf_(&m, &k, A, &m, tau, work, &neg_one, &info);
    lwork = (int) work[0];
    cgeqlf_(&m, &k, A, &m, tau, work, &neg_one, &info);
    if (lwork < (int) work[0]) lwork = (int) work[0];
    cgerqf_(&k, &m, A, &k, tau, work, &neg_one, &info);
    if (lwork < (int) work[0]) lwork = (int) work[0];
    cgelqf_(&k, &m, A, &k, tau, work, &neg_one, &info);
    if (lwork < (int) work[0]) lwork = (int) work[0];
    // free the workspace
    free(work);
    // allocate the workspace to be as large as we need
    work = (float complex *) malloc(sizeof(float complex) * lwork);
    //----------------------------------------------------------------------------------------------
    // qr
    //----------------------------------------------------------------------------------------------
    setAfc(m, k, A);
    cgeqrf_(&m, &k, A, &m, tau, work, &lwork, &info);
    //----------------------------------------------------------------------------------------------
    // AOCL
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    cungqr_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[0][0] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // Old reference
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    cungqr_ref_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[0][1] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // My version
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    cungqr_new_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[0][2] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // ql
    //----------------------------------------------------------------------------------------------
    setAfc(m, k, A);
    cgeqlf_(&m, &k, A, &m, tau, work, &lwork, &info);
    //----------------------------------------------------------------------------------------------
    // AOCL
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    cungql_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[1][0] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // Old reference
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    cungql_ref_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[1][1] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // My version
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    cungql_new_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[1][2] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // rq
    //----------------------------------------------------------------------------------------------
    setAfc(m, k, A);
    cgerqf_(&k, &m, A, &k, tau, work, &lwork, &info);
    //----------------------------------------------------------------------------------------------
    // AOCL
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    cungrq_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[2][0] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // Old reference
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    cungrq_ref_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[2][1] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // My version
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    cungrq_new_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[2][2] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // lq
    //----------------------------------------------------------------------------------------------
    setAfc(m, k, A);
    cgelqf_(&n, &m, A, &n, tau, work, &lwork, &info);
    //----------------------------------------------------------------------------------------------
    // AOCL
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    cunglq_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[3][0] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // Old reference
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    cunglq_ref_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[3][1] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // My version
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    cunglq_new_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[3][2] = elapsed_refL;
}

void timeCDouble(int m, int n, int k, double timeVals[4][3]) {
    struct timeval tp;

    double elapsed_refL;

    complex double *A = (complex double *) malloc(sizeof(complex double)*m*k);
    complex double *Q = (complex double *) malloc(sizeof(complex double)*m*n);
    complex double *tau = (complex double *) malloc(sizeof(complex double) * (m < n ? n : n)); 
    complex double *work = (complex double *) malloc(sizeof(complex double));

    int neg_one = -1;
    int lwork = -1;
    int info = 0;
    // set all of timeVals to -1 (to ensure everything is set on exit)
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            timeVals[i][j] = -1;
        }
    }
    // workspace queries for all factorization routines
    zgeqrf_(&m, &k, A, &m, tau, work, &neg_one, &info);
    lwork = (int) work[0];
    zgeqlf_(&m, &k, A, &m, tau, work, &neg_one, &info);
    if (lwork < (int) work[0]) lwork = (int) work[0];
    zgerqf_(&k, &m, A, &k, tau, work, &neg_one, &info);
    if (lwork < (int) work[0]) lwork = (int) work[0];
    zgelqf_(&k, &m, A, &k, tau, work, &neg_one, &info);
    if (lwork < (int) work[0]) lwork = (int) work[0];
    // free the workspace
    free(work);
    // allocate the workspace to be as large as we need
    work = (complex double *) malloc(sizeof(complex double) * lwork);
    //----------------------------------------------------------------------------------------------
    // qr
    //----------------------------------------------------------------------------------------------
    setAdc(m, k, A);
    zgeqrf_(&m, &k, A, &m, tau, work, &lwork, &info);
    //----------------------------------------------------------------------------------------------
    // AOCL
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    zungqr_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[0][0] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // Old reference
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    zungqr_ref_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[0][1] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // My version
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    zungqr_new_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[0][2] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // ql
    //----------------------------------------------------------------------------------------------
    setAdc(m, k, A);
    zgeqlf_(&m, &k, A, &m, tau, work, &lwork, &info);
    //----------------------------------------------------------------------------------------------
    // AOCL
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    zungql_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[1][0] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // Old reference
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    zungql_ref_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[1][1] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // My version
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    zungql_new_(&m, &n, &k, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[1][2] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // rq
    //----------------------------------------------------------------------------------------------
    setAdc(m, k, A);
    zgerqf_(&k, &m, A, &k, tau, work, &lwork, &info);
    //----------------------------------------------------------------------------------------------
    // AOCL
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    zungrq_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[2][0] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // Old reference
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    zungrq_ref_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[2][1] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // My version
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    zungrq_new_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[2][2] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // lq
    //----------------------------------------------------------------------------------------------
    setAdc(m, k, A);
    zgelqf_(&k, &m, A, &k, tau, work, &lwork, &info);
    //----------------------------------------------------------------------------------------------
    // AOCL
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    zunglq_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[3][0] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // Old reference
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    zunglq_ref_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[3][1] = elapsed_refL;
    //----------------------------------------------------------------------------------------------
    // My version
    //----------------------------------------------------------------------------------------------
    for (size_t i = 0; i < m*k; ++i) {
        Q[i] = A[i];
    }
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    zunglq_new_(&n, &m, &k, Q, &n, tau, work, &lwork, &info);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    timeVals[3][2] = elapsed_refL;
}

void printTimeArray(double timeVals[4][3]) {
    const char algNames[4][3] = {"QR", "QL", "RQ", "LQ"};
    const char impNames[3][10] = {"AOCL", "Reference", "New"};

    for (int i = 0; i < 4; ++i) {
        //printf("---------------------------------------------\n");
        for (int j = 0; j < 3; ++j) {
            printf("%s %s: %17.16e\n", algNames[i],impNames[j], timeVals[i][j]);
        }
    }
    //printf("---------------------------------------------\n");
}

void parseInputs(int *m, int *n, int *k, int argc, char *argv[]) {
    for(int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-m") == 0) {
            *m = atoi(argv[i+1]);
        } else if (strcmp(argv[i], "-n") == 0) {
            *n = atoi(argv[i+1]);
        } else if (strcmp(argv[i], "-k") == 0) {
            *k = atoi(argv[i+1]);
        }
    }
}

int main(int argc, char *argv[]){
    int m,n,k;
    m = 30;
    n = 30;
    k = 30;
    parseInputs(&m, &n, &k, argc, argv);
    double realTimeVals[4][3];
    double doubleTimeVals[4][3];

    timeReal(m,n,k,realTimeVals);
    printf("Single Precision\n");
    printTimeArray(realTimeVals);

    timeDouble(m,n,k,doubleTimeVals);
    printf("Double Precision\n");
    printTimeArray(doubleTimeVals);

    timeCReal(m,n,k,realTimeVals);
    printf("Single Complex Precision\n");
    printTimeArray(realTimeVals);

    timeCDouble(m,n,k,doubleTimeVals);
    printf("Double Complex Precision\n");
    printTimeArray(doubleTimeVals);
}
