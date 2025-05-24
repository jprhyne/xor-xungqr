#include<math.h>
#include<complex.h>
#include<stddef.h>
double myNorm(size_t m, size_t n, double complex *A) {
    double retVal = 0.0;
    for (size_t i = 0; i < m*n; ++i) {
        retVal += A[i]*conj(A[i]);
    }
    return sqrt(retVal);
}
