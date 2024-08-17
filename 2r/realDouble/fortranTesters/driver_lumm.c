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
    int info, ldl, ldu, m, n, k, lwork, nb, i, j;
    // double variables

    n = 20;
    ldl = -1;
    ldu = -1;

    for(i = 1; i < argc; ++i){
        if( strcmp( *(argv + i), "-ldu") == 0) {
            ldu  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-ldl") == 0) {
            ldl  = atoi( *(argv + i + 1) );
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

    if( ldl < 0 ) ldl = n;
    if( ldu < 0 ) ldu = n;

    char *source = SOURCE;
    printf("lumm %s:n = %4d, ldl = %4d, ldu = %4d\n",source, n, ldl, ldu);

    // Call the test file
    test_lumm_(&n, &ldl, &ldu);

    return 0;

}
