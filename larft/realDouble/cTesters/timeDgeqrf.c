/*
 * Note: this file is only for timing execution, it does NOT do 
 * any error checking. To ensure accuracy of these routines,
 * see the fortranTesters directory and relevant testers there.
 */
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
double computeDgeqrfPerf(double time, double m, double n)
{
    // Taken from LAWN 41, the cost of dlarft should not change when using 
    // the recursive definition
    double numer = 2*m*n*n - 2. / 3. *n*n*n + m*n + n*n + 14./3.*n;
    return numer / (time * 1.0e+9);
}
int main(int argc, char *argv[])
{
    // Character variables used for FORTRAN routine calls
    // integer variables for FORTRAN routine calls
    int m, n, nb, workQuery, lwork, info;
    // double variables for FORTRAN routine calls
    // matrix variables
    double *Q, *Qs, *work, *tau;
    // struct for help with timing
    struct timeval tp;
    // Local integer variables
    int i;
    // Local double variables
    double elapsed_refL;
    // parse user input.
    m = -1;
    n = -1;
    nb = -1;
    workQuery = -1;
    for(i = 1; i < argc; ++i) {
        if( strcmp( *(argv + i), "-m") == 0) {
            m  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-n") == 0) {
            n  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-nb") == 0) {
            nb  = atoi( *(argv + i + 1) );
            i++;
        }
    }
    // Make sure all inputs were provided
    if (m < 0 || n < 0 || nb < 0) {
        printf("Invalid input values. Must provide all inputs\n");
        printf("You provided: m=%d n=%d nb=%d\n",m,n,nb);
        return 1;
    }
    // make sure that m >= n
    if (m < n) {
        m = n;
    }
    // Construct our matrices and initialize them as the 0 matrix
    Q = (double *) calloc(m*n,sizeof(double));
    Qs = (double *) calloc(m*n,sizeof(double));
    if (Q == NULL) {
        printf("Failed to allocate memory for Q\n");
        return 2;
    }
    if (Qs == NULL) {
        printf("Failed to allocate memory for Qs\n");
        return 3;
    }

    // Fill up Q and Qs with the same random values.
    double tmpVal;
    for( i = 0; i < m*n; i++ ) {
        tmpVal = rand() / (double) (RAND_MAX) - 0.5e+00;
        Q[i] = tmpVal;
        Qs[i]= tmpVal;
    }
    // At this point, we have that Q = Qs
    // Create the work array to do workspace queries
    work = (double *) malloc(sizeof(double));
    // allocate the tau vector
    tau = (double *) malloc(n * sizeof(double));
    if (tau == NULL) {
        printf("Failed to allocate memory for tau\n");
        return 4;
    }
    // Determine how much workspace is needed for our operations
    dgeqrf_ref_(&m, &n, &nb, Q, &m, tau, work, &workQuery, &info );
    lwork = work[0];

    // reallocate work to be of the right size
    work = (double *) realloc(work, lwork * sizeof(double));
    if (work == NULL) {
        printf("Failed to allocate memory for work\n");
        return 5;
    }

    // Now we start the timer
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    // Call reference dgeqrf
    dgeqrf_ref_(&m, &n, &nb, Q, &m, tau, work, &lwork, &info);
    // grab the execution time
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    // Store this value 
    double refTime = elapsed_refL;
    // Compute the flop count
    double refFlop = computeDgeqrfPerf(refTime, m, n);

    // Now we start the timer
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    // Call dgeqrf with my dlarft slotted in
    dgeqrf_rec_(&m, &n, &nb, Qs, &m, tau, work, &lwork, &info);
    // grab the execution time
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    // Store this value 
    double recTime = elapsed_refL;
    // Compute the flop count
    double recFlop = computeDgeqrfPerf(recTime, m, n);

    // Print out times and flop counts
    printf("m=%d, n=%d\n", m, n);
    printf("ref:%6.4e|%6.4e\n",refTime,refFlop);
    printf("rec:%6.4e|%6.4e\n",recTime,recFlop);
}
