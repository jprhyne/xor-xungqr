/**
 * According to LAWN 41, we have the same computations for both dorgqr and dorgql
 * So the x here means either r or l
 */
double computeZungqxPerf(double m, double n, double k, double execTime) {
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
 * According to LAWN 41, we have the same computations for both dorgrq and dorglq
 * So the x here means either r or l
 */
double computeZungxqPerf(double m, double n, double k, double execTime) {
    double perfVal = 0.0;
    // Taken from lawn 41
    // numOps = 4mnk - 2(m+n)k^2 + 4/3k^3 + 2mk - k^2 - 1/3k
    // When m\neq n=k, this simplifies to
    // numOps = 2nm^2 - 2/3m^3 + m^2 -1/3m
    // When m=n\neq k, this simplifies to
    // 4m^2k - 4mk^2 + 4/3k^3 + 2mk - k^2 - 1/3k
    double numOps = 4*m*n*k - 2*(m+n)*k*k + 4./3.*k*k*k + 2*m*k - k*k - 1./3.*k;
    return numOps / (execTime * 1.0e+9);
}
