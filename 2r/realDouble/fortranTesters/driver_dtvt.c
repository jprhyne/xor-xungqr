// Testing the functionality of the tvt subroutine
// We need to ensure 
// 1) That the output of TVT is such that T = T*V**T
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
    int ldv, ldt, n, i;

    n = 20;
    ldv = -1;
    ldt = -1;

    for(i = 1; i < argc; ++i){
        if( strcmp( *(argv + i), "-ldt") == 0) {
            ldt  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-ldv") == 0) {
            ldv  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-n") == 0) {
            n  = atoi( *(argv + i + 1) );
            i++;
        }
    }

    if( ldv < 0 ) ldv = n;
    if( ldt < 0 ) ldt = n;

    char *source = SOURCE;
    printf("dtvt %s:n = %4d, ldv = %4d, ldt = %4d\n",source, n, ldv, ldt);

    // Call the fortran tester
    test_dtvt_(&n, &ldv, &ldt);
}
