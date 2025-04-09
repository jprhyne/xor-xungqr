#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
int main(int argc, char *argv[]) {
    // integer variables
    int m, n, i, j;
    // double variables
    double alphas[3], betas[3];

    // Default values
    m = 30;
    n = 20;
    alphas[0] = 0;
    betas[0] = 0;
    alphas[1] = 1;
    betas[1] = 1;
    alphas[2] = (double) rand() / (double) (RAND_MAX) - 0.5e+00;
    betas[2] = (double) rand() / (double) (RAND_MAX) - 0.5e+00;

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

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++){
            printf("dtrmmoop: alpha=%lf, beta=%lf\n", alphas[i], betas[j]);
            testdtrmmoop_(&m, &n, alphas+i, betas+j);
        }
    }

    return 0;

}
