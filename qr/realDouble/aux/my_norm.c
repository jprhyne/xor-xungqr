double my_norm(int m, int n, double *A, int lda)
{
    double ret = 0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double tmp = A[i+j*lda];
            ret += tmp * tmp;
        }
    }
    return ret;
}
