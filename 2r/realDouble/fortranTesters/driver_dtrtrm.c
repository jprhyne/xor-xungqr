#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
int main(int argc, char *argv[]) {
    // integer variables
    int n, i, j;
    // double variables
    double alphas[3];

    // Default values
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

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++){
            printf("dtrmmoop: alpha=%lf\n", alphas[i]);
            testdtrtrm_(&n, alphas+i);
        }
    }

    return 0;

}
