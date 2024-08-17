// We get exact accuracy when comparing against reference LAPACK
// but when using AOCL we get errors as large as 1e-10. Not sure why
// Probably just a different dtrmm algorithm. May investigate? 
#ifdef USE_AOCL
    #define SOURCE "AOCL"
#else
    #define SOURCE "REF_LAPACK"
#endif

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
int main(int argc, char *argv[]) {
    // integer variables
    int info, lda, ldt, m, n, k, lwork, nb, i, j;
    // double variables

    m = 30;
    n = 20;
    lda = -1;
    ldt = -1;

    for(i = 1; i < argc; ++i){
        if( strcmp( *(argv + i), "-ldt") == 0) {
            ldt  = atoi( *(argv + i + 1) );
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
    }

    if( lda < 0 ) lda = n;
    if( ldt < 0 ) ldt = m;
    if( ldt < n ) ldt = n;
    if( ldt < m ) ldt = m;
    if( lda < m ) lda = m;
    if( lda < n ) lda = n;

    char *source = SOURCE;
    printf("dtrmm vs dtrmmoop %s: m = %4d, n = %4d, lda = %4d, ldt = %4d\n",source, m, n, lda, ldt);

    // Call the test file
    testdtrmmoop_(&m,&n,&lda,&ldt);

    return 0;

}
