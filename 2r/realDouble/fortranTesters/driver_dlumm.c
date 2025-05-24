#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
int main(int argc, char *argv[]) {
    // integer variables
    int n, i;
    // double variables
    double alphas[3];
    // default values
    n = 20;
    alphas[0] = 0;
    alphas[1] = 1;
    alphas[2] = (double) rand() / (double) (RAND_MAX) - 0.5e+00;
    for(i = 1; i < argc; ++i){
        if( strcmp( *(argv + i), "-n") == 0) {
            n  = atoi( *(argv + i + 1) );
            i++;
        }
    }
    // Call the test file for each value of alpha
    for (i = 0; i < 3; ++i) {
        printf("dlumm: alpha=%lf\n", alphas[i]);
        test_dlumm_(&n, alphas+i);
    }
    return 0;

}
