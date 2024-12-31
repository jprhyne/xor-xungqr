#include<stdlib.h>
#include<math.h>
#include"myNorm.h"
// Need to include stdlib.h and myNorm
double computeQRRepresError(size_t mS, size_t nS, double *A, double *Q, double *R) {
    double alpha = 1.0;
    char aChar = 'A';
    char nChar = 'N';
    char rChar = 'R';
    char uChar = 'U';
    int m = mS;
    int n = nS;
    int dummyVal = 0;
    // Compute our representation error
    // which is given by ||A - Q*R|| / ||A||
    //
    // Compute the norm of A
    double normA;
    normA = myNorm(mS, nS, A);

    // Compute A - Q*R
    // Allocate a workspace to use
    double *work = (double *) calloc(mS*nS,sizeof(double));
    //  First: work = Q*R
    dlacpy_(&aChar, &m, &n, Q, &m, work, &m, dummyVal);
    dtrmm_(&rChar, &uChar, &nChar, &nChar, &m, &n, &alpha, R, &n, work, &m, dummyVal, dummyVal, dummyVal, dummyVal);

    // Now: work = A - work
    for (size_t i = 0; i < mS*nS; ++i) {
        work[i] = A[i] - work[i];
    }
    double represErr = myNorm(mS, nS, work);
    free(work);

    // Compute the relative error
    represErr /= normA;
    return represErr;
}

double computeRQRepresError(size_t mS, size_t nS, double *A, double *Q, double *R) {
    double alpha = 1.0;
    char aChar = 'A';
    char lChar = 'L';
    char nChar = 'N';
    char uChar = 'U';
    int m = mS;
    int n = nS;
    int dummyVal = 0;
    // Compute our representation error
    // which is given by ||A - R*Q|| / ||A||
    //
    // Compute the norm of A
    double normA;
    normA = myNorm(mS, nS, A);

    // Compute A - R*Q
    // Allocate a workspace to use
    double *work = (double *) calloc(mS*nS,sizeof(double));
    //  First: work = R*Q
    dlacpy_(&aChar, &m, &n, Q, &m, work, &m, dummyVal);
    dtrmm_(&lChar, &uChar, &nChar, &nChar, &m, &n, &alpha, R, &m, work, &m, dummyVal, dummyVal, dummyVal, dummyVal);

    // Now: work = A - work
    for (size_t i = 0; i < mS*nS; ++i) {
        work[i] = A[i] - work[i];
    }
    double represErr = myNorm(mS, nS, work);
    free(work);

    // Compute the relative error
    represErr /= normA;
    return represErr;
}

double computeQLRepresError(size_t mS, size_t nS, double *A, double *Q, double *L) {
    double alpha = 1.0;
    char aChar = 'A';
    char lChar = 'L';
    char nChar = 'N';
    char rChar = 'R';
    int m = mS;
    int n = nS;
    int dummyVal = 0;
    // Compute our representation error
    // which is given by ||A - Q*L|| / ||A||
    //
    // Compute the norm of A
    double normA;
    normA = myNorm(mS, nS, A);

    // Compute A - Q*L
    // Allocate a workspace to use
    double *work = (double *) calloc(mS*nS,sizeof(double));
    //  First: work = Q*L
    dlacpy_(&aChar, &m, &n, Q, &m, work, &m, dummyVal);
    dtrmm_(&rChar, &lChar, &nChar, &nChar, &m, &n, &alpha, L, &n, work, &m, dummyVal, dummyVal, dummyVal, dummyVal);

    // Now: work = A - work
    for (size_t i = 0; i < mS*nS; ++i) {
        work[i] = A[i] - work[i];
    }
    double represErr = myNorm(mS, nS, work);
    free(work);

    // Compute the relative error
    represErr /= normA;
    return represErr;
}

double computeLQRepresError(size_t mS, size_t nS, double *A, double *Q, double *L) {
    double alpha = 1.0;
    char aChar = 'A';
    char lChar = 'L';
    char nChar = 'N';
    int m = mS;
    int n = nS;
    int dummyVal = 0;
    // Compute our representation error
    // which is given by ||A - L*Q|| / ||A||
    //
    // Compute the norm of A
    double normA;
    normA = myNorm(mS, nS, A);

    // Compute A - L*Q
    // Allocate a workspace to use
    double *work = (double *) calloc(mS*nS,sizeof(double));
    //  First: work = L*Q
    dlacpy_(&aChar, &m, &n, Q, &m, work, &m, dummyVal);
    dtrmm_(&lChar, &lChar, &nChar, &nChar, &m, &n, &alpha, L, &m, work, &m, dummyVal, dummyVal, dummyVal, dummyVal);

    // Now: work = A - work
    for (size_t i = 0; i < mS*nS; ++i) {
        work[i] = A[i] - work[i];
    }
    double represErr = myNorm(mS, nS, work);
    free(work);

    // Compute the relative error
    represErr /= normA;
    return represErr;
}

double computeTallOrthError(size_t mS, size_t nS, double *Q) {
    int dummyVal = 0;
    // Allocate our workspace to be of size n\times  n
    double *work = (double *) calloc(nS*nS,sizeof(double));
    // Compute the orthogonality error
    double alpha= 1.0;
    double beta = 0.0;
    char tChar = 'T';
    char uChar = 'U';
    // Compute work = Q'Q - I
    int m = mS;
    int n = nS;
    dsyrk_(&uChar, &tChar, &n, &m, &alpha, Q, &m, &beta, work, &n, dummyVal, dummyVal);
    for (size_t i = 0; i < nS; ++i) 
        work[i*nS + i] -= 1.0;
    // Compute the norm of work
    double orthErr = myNorm(nS, nS, work);
    orthErr /= sqrt((double) n);
    free(work);

    return orthErr;
}

double computeWideOrthError(size_t mS, size_t nS, double *Q) {
    int dummyVal = 0;
    // Allocate our workspace to be of size m\times m (m \leq n)
    double *work = (double *) calloc(mS*mS,sizeof(double));
    // Compute the orthogonality error
    double alpha= 1.0;
    double beta = 0.0;
    char nChar = 'N';
    char uChar = 'U';
    // Compute work = QQ' - I
    int m = mS;
    int n = nS;
    dsyrk_(&uChar, &nChar, &m, &n, &alpha, Q, &m, &beta, work, &m, dummyVal, dummyVal);
    for (size_t i = 0; i < mS; ++i) 
        work[i*mS + i] -= 1.0;

    // Compute the norm of work
    double orthErr = myNorm(mS, mS, work);
    orthErr /= sqrt((double) m);
    free(work);

    return orthErr;
}
