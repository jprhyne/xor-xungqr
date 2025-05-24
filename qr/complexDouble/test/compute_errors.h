#ifndef COMPUTE_ERRORS_H
#define COMPUTE_ERRORS_H
double computeQRRepresError(size_t mS, size_t nS, double complex *A, double complex *Q, double complex *R);
double computeRQRepresError(size_t mS, size_t nS, double complex *A, double complex *Q, double complex *R);
double computeQLRepresError(size_t mS, size_t nS, double complex *A, double complex *Q, double complex *L);
double computeLQRepresError(size_t mS, size_t nS, double complex *A, double complex *Q, double complex *L);
double computeTallOrthError(size_t mS, size_t nS, double complex *Q);
double computeWideOrthError(size_t mS, size_t nS, double complex *Q);
#endif
