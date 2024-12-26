/**
 * This file is responsible for checking the accuracy of the dorgqr routine. This implicitly
 * is checking the following
 * 1) DLARFT (T computation)
 * 2) First application of H
 * 3) DLARFB0C2 (Application of H inside the loop)
 * 4) General structure of the dorgqr algorithm
 *
 * Testing in the larft directory reasonably rules out #1 from being a source of error
 * Testing in this directory of dlarfb0c2 reasonably rules out #3
 * So, if we are not accurate, it is either #2 or #4
 *
 * Note: We also compute the time of each dorgqr call, but only report it if the user supplies 
 * the '-t' flag
 */
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

void usage() {
    printf("./test_dorgqr.exe [-m numRows -n numCols -t -e -c]\n");
    printf(" -m numCols is the number of columns in  the generated matrix\n");
    printf(" -n numRows is the number of rows in  the generated matrix\n");
    printf(" -t flag that prints out just the timing information\n");
    printf(" -e flag that prints out just the error information\n");
    printf(" -c flag that prints out all data in a machine readable format\n");

}
double computeDorgqrPerf(double m, double n, double k, double execTime) {
    double perfVal = 0.0;
    // Taken from lawn 41
    // numOps = 4mnk - 2(m+n)k^2 + 4/3k^3 + 3nk - mk - k^2 - 4/3k
    // When m\neq n=k, this simplifies to
    // numOps = 2mn^2 + 1/3 n^3 + n^2 - mn - 4/3 n
    // When m=n\neq k, this simplifies to
    // numOps = 4m^2k - 4mk^2 + 4/3/ k^3 + 2mk - k^2 - 4/3 k
    double numOps = 4*m*n*k - 2*(m+n)*k*k + 4./3.*k*k*k + 3*n*k - m*k - k*k - 4./3.*k;
    return numOps / (execTime * 1.0e+9);
}

double myNorm(size_t m, size_t n, double *A) {
    double retVal = 0.0;
    for (size_t i = 0; i < m*n; ++i) {
        retVal += A[i]*A[i];
    }
    return sqrt(retVal);
}

int main(int argc, char *argv[]) {
    // Scalars
    // integers
    int info, ldt, m, n, lwork, i, j, neg_one;
    // size_t (m*n*sizeof(double) may be larger than 
    // can be represented in an int, but these fortran subroutines need signed 
    // integers, so we keep these as separate for allocating data
    size_t mS, nS, lworkS;
    // double
    double norm_orth, norm_repres, elapsed, alpha, beta;
    // logical
    bool printTimes, printErrors, printCompact;
    // character
    char fChar, nChar, rChar, uChar;
    // Arrays
    // double
    double *A, *As, *R, *tau, *work;
    // Structs
    // timing helper
    struct timeval tp;

    // Set default values
    printTimes = false;
    printErrors = false;
    printCompact = false;
    m = 30;
    n = 20;
    fChar = 'F';
    nChar = 'N';
    rChar = 'R';
    uChar = 'U';
    // dummy value for use with functions that have character inputs
    int dummyVal = 0;

    // Parse input flags
    for(i = 1; i < argc; ++i) {
        // Check for the number of rows
        if( strcmp( argv[i], "-m" ) == 0) {
            m = atoi(argv[i+1]);
            i++;
        }
        // Check for the number of columns
        if( strcmp( argv[i], "-n" ) == 0) {
            n = atoi(argv[i+1]);
            i++;
        }
        // Check if the user wants to print out timing information or not
        if( strcmp( argv[i], "-t" ) == 0) {
            printTimes = true;
        }

        if( strcmp( argv[i], "-e" ) == 0) {
            printErrors = true;
        }

        if( strcmp( argv[i], "-c" ) == 0) {
            printCompact = true;
        }
    }

    // sets size_t variables to be the same value as their integer counterparts
    mS = (size_t) m;
    nS = (size_t) n;
    // Allocate our arrays
    A  = (double *) malloc(mS*nS*sizeof(double));
    As = (double *) malloc(mS*nS*sizeof(double));
    tau  = (double *) malloc(mS*sizeof(double));
    work = (double *) malloc(sizeof(double));
    // Allocate a helper array R to be all 0's
    R = (double *) calloc(nS*nS,sizeof(double));

    neg_one = -1;
    // Fill our arrays with random values (and copy A into As at the same time)
    double tmpVal;
    for( i = 0; i < mS*nS; ++i ) {
        tmpVal = rand() / (double) (RAND_MAX) - 0.5e+00;
        A[i] = tmpVal;
        As[i]= tmpVal;
    }
    for ( i = 0; i < mS; ++i ) {
        tmpVal = rand() / (double) (RAND_MAX) - 0.5e+00;
        tau[i] = tmpVal;
    }
    // Now, we start the process.
    // Perform our workspace query
    // We will need to factorize then construct Q. So we must query
    // dgeqrf and dorgqr
    dgeqrf_(&m, &n, A, &m, tau, work, &neg_one, &info);
    // Store this value
    lwork = (int) work[0];
    // Query reference dorgqr as this is the baseline and workspace requirements
    // will only lower from here not increase.
    dorgqr_ref_(&m, &n, &n, A, &m, tau, work, &neg_one, &info);
    
    if (lwork < (int) work[0]) lwork = (int) work[0];

    free(work);
    lworkS = (size_t) lwork;
    work = (double *) malloc(lworkS * sizeof(double));

    // Factorize A into [R]
    //                  [Q]
    dgeqrf_(&m, &n, A, &m, tau, work, &lwork, &info);

    // Store the R component of A in R
    dlacpy_(&uChar, &m, &n, A, &m, R, &n, dummyVal);

    // Start the timer 
    gettimeofday(&tp,NULL);
    elapsed=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    // Construct Q
    dorgqr_ref_(&m, &n, &n, A, &m, tau, work, &lwork, &info);

    // Stop the timer
    gettimeofday(&tp,NULL);
    elapsed+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    // Store the reference value
    double refTime = elapsed;

    // Compute the performance metric.
    double refPerf = computeDorgqrPerf(m, n, n, refTime);

    // Now, we compute the accuracy metrics. The first is representation error
    // which is given by ||As - Q*R|| / ||As||
    //
    // Compute the norm of A
    double normA;
    normA = myNorm(mS, nS, As);

    // Compute A - Q*R
    //
    // Allocate a workspace to use
    free(work);
    work = (double *) malloc(mS*nS*sizeof(double));
    //  First: work = Q*R
    dlacpy_(&fChar, &m, &n, A, &m, work, &m, NULL, dummyVal);
    alpha = 1.0;
    dtrmm_(&rChar, &uChar, &nChar, &nChar, &m, &n, &alpha, R, &n, work, &m, dummyVal, dummyVal, dummyVal, dummyVal);

    // Now: work = A - work
    tmpVal = 0.0;
    for (i = 0; i < mS*nS; ++i) {
        work[i] = As[i] - work[i];
        tmpVal += work[i] * work[i];
    }

    // Compute the norm of work
    double relativeErr = sqrt(tmpVal);
    // Compute the relative error
    relativeErr /= normA;

    // Print this information to the console along with some diagnostics
    if (printCompact) {
        // Say what the printed output will look like
        if (!printTimes && !printErrors) 
            printf("Only sources will be printed");
        else 
            printf("source");
        if (printTimes)
            printf(":perf");
        if (printErrors)
            printf(":error");
        printf("\n");
    }
    if (printCompact) {
        printf("ref");
    } else {
        printf("m=%d, n=%d, printTimings=%d, printErrors=%d\n", m, n, printTimes, printErrors);
    }
    if (printTimes) {
        if( printCompact ) {
            printf(":%6.4e", refPerf);
        } else {
            printf("reference performance: %6.4e\n", refPerf);
        }
    }
    if (printErrors) {
        if( printCompact ) {
            printf(":%6.4e", relativeErr);
        } else {
            printf("reference forward error: %6.4e\n", relativeErr);
        }
    }
    if(printCompact){
        printf("\n");
    }

    // free our arrays
freeMemory:
    free(A);
    free(As);
    free(tau);
    free(work);
    free(R);
}
