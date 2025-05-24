// This header file is to declare the mangled fortran lapack functions for use in test_zungz.c
#ifndef MY_LAPACK_HEADERS_H
#define MY_LAPACK_HEADERS_H
    // QR
    void zgeqrf_(int *m, int *n, double complex *A, int *lda, double complex *tau, double complex *work, int *lwork, int *info);
    void zungqr_ref_(int *m, int *n, int *k, double complex *A, int *lda, double complex *tau, double complex *work, int *lwork, int *info);
    void zungqr_(int *m, int *n, int *k, double complex *A, int *lda, double complex *tau, double complex *work, int *lwork, int *info);
    void zungqr_zlarfb0c2_(int *m, int *n, int *k, double complex *A, int *lda, double complex *tau, double complex *work, int *lwork, int *info);
    // LQ
    void zgelqf_(int *m, int *n, double complex *A, int *lda, double complex *tau, double complex *work, int *lwork, int *info);
    void zunglq_ref_(int *m, int *n, int *k, double complex *A, int *lda, double complex *tau, double complex *work, int *lwork, int *info);
    void zunglq_(int *m, int *n, int *k, double complex *A, int *lda, double complex *tau, double complex *work, int *lwork, int *info);
    void zunglq_zlarfb0c2_(int *m, int *n, int *k, double complex *A, int *lda, double complex *tau, double complex *work, int *lwork, int *info);
    // QL
    void zgeqlf_(int *m, int *n, double complex *A, int *lda, double complex *tau, double complex *work, int *lwork, int *info);
    void zungql_ref_(int *m, int *n, int *k, double complex *A, int *lda, double complex *tau, double complex *work, int *lwork, int *info);
    void zungql_(int *m, int *n, int *k, double complex *A, int *lda, double complex *tau, double complex *work, int *lwork, int *info);
    void zungql_zlarfb0c2_(int *m, int *n, int *k, double complex *A, int *lda, double complex *tau, double complex *work, int *lwork, int *info);
    // RQ
    void zgerqf_(int *m, int *n, double complex *A, int *lda, double complex *tau, double complex *work, int *lwork, int *info);
    void zungrq_ref_(int *m, int *n, int *k, double complex *A, int *lda, double complex *tau, double complex *work, int *lwork, int *info);
    void zungrq_(int *m, int *n, int *k, double complex *A, int *lda, double complex *tau, double complex *work, int *lwork, int *info);
    void zungrq_zlarfb0c2_(int *m, int *n, int *k, double complex *A, int *lda, double complex *tau, double complex *work, int *lwork, int *info);
#endif
