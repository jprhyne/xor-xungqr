#include<math.h>
#include<stddef.h>
double myNorm(size_t m, size_t n, double *A) {
    double retVal = 0.0;
    for (size_t i = 0; i < m*n; ++i) {
        retVal += A[i]*A[i];
    }
    return sqrt(retVal);
}
