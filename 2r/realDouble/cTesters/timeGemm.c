#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
// This function times calls to gemm in order to 
// give a metric for what kind of performance is expected
int main(int argc, char *argv[]) 
{
    long n = 100;
    double *A, *B, *C = NULL;
    double perform_refL, elapsed_refL;
    char nChar = 'N';
    double one =1.0;
    double zero=0.0;
    size_t dummy = 0;
    // struct for help with timing
    struct timeval tp;
    for(int i = 1; i < argc; ++i){
        if( strcmp( *(argv + i), "-n") == 0) {
            n  = atoi( *(argv + i + 1) );
            i++;
        }
    }
    A = (double *) malloc(n*n*sizeof(double));
    B = (double *) malloc(n*n*sizeof(double));
    C = (double *) malloc(n*n*sizeof(double));
    for (long i = 0; i < n * n; i++) {
        A[i] = (double) rand() / (double) (RAND_MAX) - 0.5e+00;
        B[i] = (double) rand() / (double) (RAND_MAX) - 0.5e+00;
        C[i] = (double) rand() / (double) (RAND_MAX) - 0.5e+00;
    }
    
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    dgemm_(&nChar, &nChar, &n, &n, &n, &one, A, &n, B, &n, &one, C, &n, dummy, dummy);
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    perform_refL = 2.0 * (double) n * (double) n * (double) n / (elapsed_refL*1.0e+9);
    printf("%2.10f\n", perform_refL);
}
