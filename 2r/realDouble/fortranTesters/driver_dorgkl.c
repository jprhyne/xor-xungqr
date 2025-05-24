#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
int main(int argc, char *argv[]) {
    // integer variables
    int info, m, n, i;
    // double variables

    m = 30;
    n = 20;

    for(i = 1; i < argc; ++i){
        if( strcmp( *(argv + i), "-m") == 0) {
            m  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-n") == 0) {
            n  = atoi( *(argv + i + 1) );
            i++;
        }
    }

    printf("dgeqrf dlarft dorgkl: m = %4d, n = %4d\n", m, n);

    // Call the test file
    test_dorgkl_(&m, &n);

    return 0;

}
