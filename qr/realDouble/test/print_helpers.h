#ifndef PRINTHELPERS_H
#define PRINTHELPERS_H
void printHeader(bool *printFlags, int m, int n);
void printInfo(double perfVal, double repErr, double orthErr, bool *printFlags, char *longSource,
        char *shortSource);
#endif
