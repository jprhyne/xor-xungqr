#ifndef COMPUTE_ERRORS_H
#define COMPUTE_ERRORS_H
double computeQRRepresError(size_t mS, size_t nS, double *A, double *Q, double *R);
double computeRQRepresError(size_t mS, size_t nS, double *A, double *Q, double *R);
double computeQLRepresError(size_t mS, size_t nS, double *A, double *Q, double *L);
double computeLQRepresError(size_t mS, size_t nS, double *A, double *Q, double *L);
double computeTallOrthError(size_t mS, size_t nS, double *Q);
double computeWideOrthError(size_t mS, size_t nS, double *Q);
#endif
