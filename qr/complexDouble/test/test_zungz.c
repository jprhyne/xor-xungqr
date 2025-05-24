/**
 * This file is responsible for checking the accuracy of the zung{qr,lq,ql,rq} routines. 
 * This implicitly is checking the following
 * 1) DLARFT (T computation)
 * 2) First application of H
 * 3) DLARFB0C2 (Application of H inside the loop)
 * 4) General structure of the zungqr algorithm
 *
 * Testing in the larft directory reasonably rules out #1 from being a source of error
 * However, any other error is being checked here. This means if you get an error, it could be
 * any of these 3 steps. 
 */
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <sys/time.h>
// myNorm
#include "myNorm.h" 
//computeQRRepresError
//computeRQRepresError
//computeQLRepresError
//computeLQRepresError
//computeTallOrthError
//computeWideOrthError
#include "compute_errors.h"
//computeZungqxPerf
//computeZungxqPerf
#include "perf_calc.h"
//printHeader
//printInfo
#include "print_helpers.h"
//zgeqrf_
//zungqr_ref_
#include "lapackHeaders.h"

void usage() {
    printf("./test_zungqr.exe [-m numRows -n numCols -t -p -e -c -v -qr -lq -ql -rq -h]\n");
    printf(" -m  numCols is the number of columns in the generated matrix\n");
    printf("\tDefault is 30\n");
    printf(" -n  numRows is the number of rows in the generated matrix\n");
    printf("\tDefault is 20\n");
    printf(" -t  flag that prints out the timing information\n");
    printf(" -p  flag that prints out the performance information\n");
    printf(" -e  flag that prints out the error information\n");
    printf(" -c  flag that prints out all data in a machine readable format\n");
    printf(" -v  flag that prints out a header at the top for compact usage\n");
    printf("\tNote: only respected if -c is given\n");
    printf(" -qr flag that states to print out qr information\n");
    printf(" -lq flag that states to print out lq information\n");
    printf(" -ql flag that states to print out ql information\n");
    printf(" -rq flag that states to print out rq information\n");
    printf("\tNote: if none of the qr, lq, ql, rq flags are given, then all are printed\n");
    printf(" -h Print this usage information and exit\n");
}

// typedef to allow passing in different factorization functions to our tester (zge{qr,lq,ql,rq}f)
typedef void (* factorizeFunc)(int *m, int *n, double complex *A, int *lda, double complex *tau, double complex *work, 
        int *lwork, int *info);

// typedef to allow easier passing in of different functions to construct Q (zung{qr,lq,ql,rq})
typedef void (* constructQFunc)(int *m, int *n, int *k, double complex *A, int *lda, double complex *tau, double complex *work,
        int *lwork, int *info);

typedef double (*computePerfFunc)(double m, double n, double k, double timeSpent);

typedef double (*computeOrthFunc)(size_t mS, size_t nS, double complex *Q);

// X represents either R or L depending on which factorization is used
typedef double (*computeRepresFunc)(size_t mS, size_t nS, double complex *A, double complex *Q, double complex *X);

// This function will compute the orthogonality and representation error metrics as well as the 
// performance metrics. The inputs are typedef'd function pointers. They are defined directly before
// this function definition.
/**
 * Inputs
 * factor:              pointer to a function that will factorize A into A = QX or A=XQ where 
 *                          X is either upper or lower triangular
 * formQ :              pointer to a function that will convert the output of factor into a Q matrix
 * computePerf:         pointer to a function that will compute the performance metric for forming Q
 * computeRepr:         pointer to a function that will compute the representation error for A = QX or A = XQ
 * mS:                  The number of rows in A
 * nS:                  The number of columns in A
 * upperFactor:         True iff our triangular factor is upper triangular [qr and rq]
 * factorStartTop:      True iff our triangular factor starts at Q(1,1) [qr and lq]
 * timeFactorization:   True iff we want to time the factorization step. This will be false almost 
 *                          all of the time
 * Outputs
 *  A pointer to a heap allocated struct containing our results. Caller is responsible for freeing
 *  this memory.
 *      retVal[0] is the timing metrics
 *      retVal[1] is the performance value
 *      retVal[2] is the representation error
 *      retVal[3] is the orthogonality error
 */
double* computeMetrics(factorizeFunc factor, constructQFunc formQ, computePerfFunc computePerf,
        computeRepresFunc computeRepr, computeOrthFunc computeOrthErr, size_t mS, size_t nS,
        double complex *A, bool upperFactor, bool factorStartTop, bool timeFactorization) {
    // Scalars
    int neg_one, m, n, lwork, info, i, j, minDim, maxDim;
    size_t lworkS, minDimS, maxDimS;
    // double complex
    double complex alpha, beta;
    // double real
    double norm_orth, norm_repres, elapsed;
    // character
    char aChar, lChar, uChar;
    // Arrays
    // double complex
    double complex *Q, *X, *tau, *work; // X is either L or R depending on the routine used.
    // return array
    double *retVal;
    // Structs
    // timing helper
    struct timeval tp;
    // Beginning of executable statements
    neg_one = -1;
    m = (int) mS;
    n = (int) nS;
    aChar = 'A';
    lChar = 'L';
    uChar = 'U';
    minDimS = (mS < nS) ? mS : nS;
    maxDimS = (minDimS == mS) ? nS : mS;
    minDim = (int) minDimS;
    maxDim = (int) maxDimS;
    elapsed = 0.0;
    // Allocate our memory
    Q = (double complex *) malloc(mS*nS*sizeof(double complex));
    X = (double complex *) malloc(minDim*minDim*sizeof(double complex));
    tau = (double complex *) malloc(minDim*sizeof(double complex));
    work = (double complex *) malloc(sizeof(double complex));
    retVal = (double *) calloc(4, sizeof(double complex));

    // Copy A into Q
    zlacpy_(&aChar, &m, &n, A, &m, Q, &m);

    // Query the factorize and formQ functions for workspace needs
    factor(&m, &n, Q, &m, tau, work, &neg_one, &info);
    lwork = ((int) work[0]);
    formQ(&m, &n, &minDim, Q, &m, tau, work, &neg_one, &info);
    if (lwork < (int) work[0]) lwork = (int) work[0];
    lworkS = (size_t) lwork;
    // reallocate our workspace to be the proper size
    free(work);
    work = (double complex *) malloc(lworkS*sizeof(double complex));
    // factorize our matrix Q into V, X
    if (timeFactorization) {
        gettimeofday(&tp,NULL);
        elapsed=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    }
    factor(&m, &n, Q, &m, tau, work, &lwork, &info);
    if (timeFactorization) {
        gettimeofday(&tp,NULL);
        elapsed+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec); 
    }
    // Copy over X from V depending on the proper input flags
    char triangle = (upperFactor) ? 'U' : 'L';
    if (factorStartTop) {
        zlacpy(&triangle, &minDim, &minDim, Q, &m, X, &minDim);
    } else {
        if (upperFactor) {
            //qToR.f
            qtor_(&m, &n, Q, X);
        } else {
            //qToL.f
            qtol_(&m, &n, Q, X);
        }
    }
    // now form Q
    gettimeofday(&tp,NULL);
    // May be timing the factorization or not here.
    elapsed+=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
    formQ(&m, &n, &minDim, Q, &m, tau, work, &lwork, &info);
    gettimeofday(&tp,NULL);
    elapsed+=((double complex)tp.tv_sec+(1.e-6)*tp.tv_usec); 
    // Now compute our metrics
    retVal[0] = elapsed;
    retVal[1] = computePerf((double) m, (double) n, (double) minDim, elapsed);
    retVal[2] = computeRepr(mS, nS, A, Q, X);
    retVal[3] = computeOrthErr(mS, nS, Q);

    // Free all our arrays
    free(Q);
    free(X);
    free(tau);
    free(work);
    // Return our metrics. Caller MUST free this (32 bytes is a lot of data you know!)
    return retVal;
}

int main(int argc, char *argv[]) {
    // Scalars
    // size_t (these are used to allocate data and thus should not be signed because -1 bytes is
    // nonsense)
    size_t m, n;
    // double complex
    double norm_orth, norm_repres, elapsed;
    // logical
    bool printTime, printPerf, printErrors, printCompact, verboseHeader;
    bool printQR, printLQ, printQL, printRQ;
    // character
    char aChar, uChar;
    // Arrays
    // double complex
    double complex *A, *At;
    // Structs

    // Set default values
    printTime = false;
    printPerf = false;
    printErrors = false;
    printCompact = false;
    verboseHeader = false;
    printQR = false;
    printLQ = false;
    printQL = false;
    printRQ = false;
    m = 30;
    n = 20;
    aChar = 'A';
    uChar = 'U';
    // dummy value for use with functions that have character inputs
    int dummyVal = 0;

    // Parse input flags
    for(int i = 1; i < argc; ++i) {
        // Check for the number of rows
        if( strcmp( argv[i], "-m" ) == 0) {
            m = (size_t) atoi(argv[i+1]);
            i++;
        }
        // Check for the number of columns
        if( strcmp( argv[i], "-n" ) == 0) {
            n = (size_t) atoi(argv[i+1]);
            i++;
        }
        // Check if the user wants to print out timing information or not
        if( strcmp( argv[i], "-t" ) == 0) {
            printTime = true;
        }

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

        if( strcmp( argv[i], "-qr" ) == 0) {
            printQR = true;
        }

        if( strcmp( argv[i], "-lq" ) == 0) {
            printLQ = true;
        }

        if( strcmp( argv[i], "-ql" ) == 0) {
            printQL = true;
        }

        if( strcmp( argv[i], "-rq" ) == 0) {
            printRQ = true;
        }

        if( strcmp( argv[i], "-h" ) == 0) {
            usage();
            return 0;
        }
    }
    // Store our printing flags in an array for ease of passing into functions
    bool flagVec[5] = {printTime, printPerf, printErrors, printCompact, verboseHeader};
    bool printAny = printTime || printPerf || printErrors;
    if (!(printQR || printLQ || printQL || printRQ)) {
        printQR = true;
        printLQ = true;
        printQL = true;
        printRQ = true;
    }
    // sets size_t variables to be the same value as their integer counterparts
    // Allocate our arrays
    A  = (double complex *) malloc(m*n*sizeof(double complex));

    // Fill our arrays with random values (and copy A into As at the same time)
    double complex tmpVal;
    for( size_t i = 0; i < m*n; ++i ) {
        tmpVal = (rand() + rand()*I) / (double complex) (RAND_MAX) - 0.5e+00;
        A[i] = tmpVal;
    }
    double *metrics = NULL;
    printHeader(flagVec, m, n);
    if(printQR) {
        //----------------------------------------------------------------------------------------//
        // Try QR                                                                                 //
        //----------------------------------------------------------------------------------------//
        //----------------------------------------------------------------------------------------//
        // Reference                                                                              //
        //----------------------------------------------------------------------------------------//
        metrics = computeMetrics(zgeqrf_, zungqr_ref_, computeZungqxPerf, 
                computeQRRepresError, computeTallOrthError, m, n, A, true /*upperFactor*/,
                true /*factorStartTop*/, false /*timeFactorization*/);
        if(printAny) printf("QR metrics\n");
        printInfo(metrics, flagVec, "reference", "ref");
        free(metrics);
        metrics = NULL;
        //----------------------------------------------------------------------------------------//
        // Optimized                                                                              //
        //----------------------------------------------------------------------------------------//
        metrics = computeMetrics(zgeqrf_, zungqr_, computeZungqxPerf, computeQRRepresError, 
                computeTallOrthError, m, n, A, true, true, false);
        printInfo(metrics, flagVec, "optimized", "opt");
        free(metrics);
        metrics = NULL;
        //----------------------------------------------------------------------------------------//
        // DLARFB0C2                                                                              //
        //----------------------------------------------------------------------------------------//
        metrics = computeMetrics(zgeqrf_, zungqr_zlarfb0c2_, computeZungqxPerf, computeQRRepresError, 
                computeTallOrthError, m, n, A, true, true, false);
        printInfo(metrics, flagVec, "DLARFB0C2", "0c2");
        free(metrics);
        metrics = NULL;
    }
    if (printLQ) {
        //----------------------------------------------------------------------------------------//
        // Try LQ                                                                                 //
        //----------------------------------------------------------------------------------------//
        //----------------------------------------------------------------------------------------//
        // Reference                                                                              //
        //----------------------------------------------------------------------------------------//
        metrics = computeMetrics(zgelqf_, zunglq_ref_, computeZungxqPerf, 
                computeLQRepresError, computeWideOrthError, n, m, A, false /*upperFactor*/,
                true /*factorStartTop*/, false /*timeFactorization*/);
        if(printAny) printf("LQ metrics\n");
        printInfo(metrics, flagVec, "reference", "ref");
        free(metrics);
        metrics=NULL;
        //----------------------------------------------------------------------------------------//
        // Optimized                                                                              //
        //----------------------------------------------------------------------------------------//
        metrics = computeMetrics(zgelqf_, zunglq_, computeZungxqPerf, 
                computeLQRepresError, computeWideOrthError, n, m, A, false /*upperFactor*/,
                true /*factorStartTop*/, false /*timeFactorization*/);
        printInfo(metrics, flagVec, "optimized", "opt");
        free(metrics);
        metrics=NULL;
        //----------------------------------------------------------------------------------------//
        // DLARFB0C2                                                                              //
        //----------------------------------------------------------------------------------------//
        metrics = computeMetrics(zgelqf_, zunglq_zlarfb0c2_, computeZungxqPerf, 
                computeLQRepresError, computeWideOrthError, n, m, A, false /*upperFactor*/,
                true /*factorStartTop*/, false /*timeFactorization*/);
        printInfo(metrics, flagVec, "DLARFB0C2", "0c2");
        free(metrics);
        metrics=NULL;
    }
    if (printQL) {
        //----------------------------------------------------------------------------------------//
        // Try QL                                                                                 //
        //----------------------------------------------------------------------------------------//
        //----------------------------------------------------------------------------------------//
        // Reference                                                                              //
        //----------------------------------------------------------------------------------------//
        metrics = computeMetrics(zgeqlf_, zungql_ref_, computeZungqxPerf, 
                computeQLRepresError, computeTallOrthError, m, n, A, false /*upperFactor*/,
                false /*factorStartTop*/, false /*timeFactorization*/);
        if(printAny) printf("QL metrics\n");
        printInfo(metrics, flagVec, "reference", "ref");
        free(metrics);
        metrics=NULL;
        //----------------------------------------------------------------------------------------//
        // Optimized                                                                              //
        //----------------------------------------------------------------------------------------//
        metrics = computeMetrics(zgeqlf_, zungql_, computeZungqxPerf, 
                computeQLRepresError, computeTallOrthError, m, n, A, false /*upperFactor*/,
                false /*factorStartTop*/, false /*timeFactorization*/);
        printInfo(metrics, flagVec, "optimized", "opt");
        free(metrics);
        metrics=NULL;
        //----------------------------------------------------------------------------------------//
        // DLARFB0C2                                                                              //
        //----------------------------------------------------------------------------------------//
        metrics = computeMetrics(zgeqlf_, zungql_zlarfb0c2_, computeZungqxPerf, 
                computeQLRepresError, computeTallOrthError, m, n, A, false /*upperFactor*/,
                false /*factorStartTop*/, false /*timeFactorization*/);
        printInfo(metrics, flagVec, "DLARFB0C2", "0c2");
        free(metrics);
        metrics=NULL;
    }
    if (printRQ) {
        //----------------------------------------------------------------------------------------//
        // Try RQ                                                                                 //
        //----------------------------------------------------------------------------------------//
        //----------------------------------------------------------------------------------------//
        // Reference                                                                              //
        //----------------------------------------------------------------------------------------//
        metrics = computeMetrics(zgerqf_, zungrq_ref_, computeZungxqPerf, 
                computeRQRepresError, computeWideOrthError, n, m, A, true /*upperFactor*/,
                false /*factorStartTop*/, false /*timeFactorization*/);
        if(printAny) printf("RQ metrics\n");
        printInfo(metrics, flagVec, "reference", "ref");
        free(metrics);
        metrics=NULL;
        //----------------------------------------------------------------------------------------//
        // Optimized                                                                              //
        //----------------------------------------------------------------------------------------//
        metrics = computeMetrics(zgerqf_, zungrq_, computeZungxqPerf, 
                computeRQRepresError, computeWideOrthError, n, m, A, true /*upperFactor*/,
                false /*factorStartTop*/, false /*timeFactorization*/);
        printInfo(metrics, flagVec, "optimized", "opt");
        free(metrics);
        metrics=NULL;
        //----------------------------------------------------------------------------------------//
        // DLARFB0C2                                                                              //
        //----------------------------------------------------------------------------------------//
        metrics = computeMetrics(zgerqf_, zungrq_zlarfb0c2_, computeZungxqPerf, 
                computeRQRepresError, computeWideOrthError, n, m, A, true /*upperFactor*/,
                false /*factorStartTop*/, false /*timeFactorization*/);
        printInfo(metrics, flagVec, "DLARFB0C2", "0c2");
        free(metrics);
        metrics=NULL;
    }
freeMemory:
    // free our arrays
    free(A);
}
