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
 */
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

void usage() {
    printf("./test_dorgqr.exe [-m numRows -n numCols -p -e -c -h]\n");
    printf(" -m numCols is the number of columns in  the generated matrix\n");
    printf(" -n numRows is the number of rows in  the generated matrix\n");
    printf(" -p flag that prints out the performance information\n");
    printf(" -e flag that prints out the error information\n");
    printf(" -c flag that prints out all data in a machine readable format\n");
    printf(" -v flag that prints out a header at the top for compact usage.\n\tNote: only respected if -c is given\n");
    printf(" -h Print this usage information and exit\n");
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
/**
 * This function prints the header of our output. This is mainly too keep the
 * main function more streamlined and easier to follow
 * printFlags: Array of length 3 that contains the following elements
 *  printFlags[0] == true means we print out the performance value
 *  printFlags[1] == true means we print out relative error
 *  printFlags[2] == true iff we are printing out in a compact (machine readable) format
 *  printFlags[3] == true iff we want to print out the header when we are in compact printing mode
 *
 * m: The number of rows in the matrix A
 * n: The number of columns in the matrix A
 */
void printHeader(bool *printFlags, int m, int n) {
    // Grab the flags
    bool printPerf    = printFlags[0];
    bool printErrors  = printFlags[1];
    bool printCompact = printFlags[2];
    bool verbosePrint = printFlags[3];
    // Print out some diagnostics for the user
    if (printCompact && verbosePrint) {
        // Say what the printed output will look like
        if (!printPerf && !printErrors) 
            printf("Only sources will be printed");
        else 
            printf("source");
        if (printPerf)
            printf(":perf");
        if (printErrors)
            printf(":repErr:orthErr");
        printf("\n");
    }
    if (!printCompact) {
        printf("m=%d, n=%d, printErrors=%d, printPerformance=%d\n", m, n, printErrors, printPerf);
    }
}
/**
 * relErr: The relative error we are printing
 * perfVal:The performance value to print
 * printFlags: Array of length 3 that contains the following elements
 *  printFlags[0] == true means we print out the performance value
 *  printFlags[1] == true means we print out relative error
 *  printFlags[2] == true iff we are printing out in a compact (machine readable) format
 *  printFlags[3] Not used in this function
 * source: string that displays the source of these values. Example
 *  values include "reference", "optimized", etc.
 */
void printInfo(double perfVal, double repErr, double orthErr, bool *printFlags, char *source) {
    // Grab the flags
    bool printPerf    = printFlags[0];
    bool printErrors  = printFlags[1];
    bool printCompact = printFlags[2];
    // Print our data in the order implied above
    // Print out the beginning of the line if we are compact. Note that 
    // if we are compact, we truncate the source string to be 3 characters long
    if (printCompact) {
        printf("%.3s",source);
    }
    // Times first
    if (printPerf) {
        if( printCompact ) {
            printf(":%6.4e", perfVal);
        } else {
            // %s is scary!
            printf("%s performance: %6.4e\n", source, perfVal);
        }
    }
    // Errors next
    if (printErrors) {
        if( printCompact ) {
            printf(":%6.4e:%6.4e", repErr, orthErr);
        } else {
            // %s is scary!
            printf("%s representation error: %6.4e. Orthogonality Error: %6.4e\n", source, repErr, orthErr);
        }
    }
    if(printCompact){
        printf("\n");
    }
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
    bool printPerf, printErrors, printCompact, verboseHeader;
    // character
    char fChar, nChar, rChar, tChar, uChar;
    // Arrays
    // double
    double *A, *As, *R, *tau, *work;
    // Structs
    // timing helper
    struct timeval tp;

    // Set default values
    printPerf = false;
    printErrors = false;
    printCompact = false;
    verboseHeader = false;
    m = 30;
    n = 20;
    fChar = 'F';
    nChar = 'N';
    rChar = 'R';
    tChar = 'T';
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
        if( strcmp( argv[i], "-p" ) == 0) {
            printPerf = true;
        }

        if( strcmp( argv[i], "-e" ) == 0) {
            printErrors = true;
        }

        if( strcmp( argv[i], "-c" ) == 0) {
            printCompact = true;
        }

        if( strcmp( argv[i], "-v" ) == 0) {
            verboseHeader = true;
        }

        if( strcmp( argv[i], "-h" ) == 0) {
            usage();
            return 0;
        }
    }
    // Store our printing flags in an array for ease of passing into functions
    bool flagVec[4] = {printPerf, printErrors, printCompact, verboseHeader};

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
    double represErr = sqrt(tmpVal);
    // Compute the relative error
    represErr /= normA;

    // Set work to 0
    for (i = 0; i < mS*nS; ++i) {
        work[i] = 0.0;
    }

    // Compute the orthogonality error
    beta = 0.0;
    // Compute work = Q'Q - I
    dsyrk_(&uChar, &tChar, &n, &m, &alpha, A, &m, &beta, work, &m);
    for (i = 0; i < nS; ++i) 
        work[i*mS + i] -= 1.0;
    // Compute the norm of work
    tmpVal = 0;
    for (i = 0; i < mS * nS; ++i)
        tmpVal += work[i] * work[i];

    double orthErr = sqrt(tmpVal);
    orthErr /= sqrt((double) n);

    // Print this information to the console along with some diagnostics
    printHeader(flagVec, m, n);
    printInfo(refPerf, represErr, orthErr, flagVec, "reference");

    // free our arrays
freeMemory:
    free(A);
    free(As);
    free(tau);
    free(work);
    free(R);
}
