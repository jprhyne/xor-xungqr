#include<stdio.h>
#include<stdbool.h>
/**
 * This function prints the header of our output. This is mainly too keep the
 * main function more streamlined and easier to follow
 * printFlags: Array of length 4 that contains the following elements
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
 * printFlags: Array of length 4 that contains the following elements
 *  printFlags[0] == true means we print out the performance value
 *  printFlags[1] == true means we print out relative error
 *  printFlags[2] == true iff we are printing out in a compact (machine readable) format
 *  printFlags[3] Not used in this function
 * longSource: string that displays the source of these values to be printed with
 * non compact printing
 *  values include "reference", "optimized", etc.
 * shortSource: string that displays the source of these values to be printed with
 * compact printing.
 *  values include "ref", "opt", "rec_T"
 */
void printInfo(double perfVal, double repErr, double orthErr, bool *printFlags, 
        char *longSource, char *shortSource) {
    // Grab the flags
    bool printPerf    = printFlags[0];
    bool printErrors  = printFlags[1];
    bool printCompact = printFlags[2];
    // Print our data in the order implied above
    // Print out the beginning of the line if we are compact.
    if (printCompact) {
        printf("%s",shortSource);
    }
    // Times first
    if (printPerf) {
        if( printCompact ) {
            printf(":%6.4e", perfVal);
        } else {
            printf("%s performance: %6.4e\n", longSource, perfVal);
        }
    }
    // Errors next
    if (printErrors) {
        if( printCompact ) {
            printf(":%6.4e:%6.4e", repErr, orthErr);
        } else {
            printf("%s representation error: %6.4e. Orthogonality Error: %6.4e\n",
                    longSource, repErr, orthErr);
        }
    }
    if(printCompact){
        printf("\n");
    }
}
