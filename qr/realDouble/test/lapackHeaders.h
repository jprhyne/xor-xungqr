// This header file is to declare the mangled fortran lapack functions for use in test_dorgz.c
#ifndef MY_LAPACK_HEADERS_H
#define MY_LAPACK_HEADERS_H
    // QR
    void dgeqrf_(int *m, int *n, double *A, int *lda, double *tau, double *work, int *lwork, int *info);
    void dorgqr_ref_(int *m, int *n, int *k, double *A, int *lda, double *tau, double *work, int *lwork, int *info);
    void dorgqr_(int *m, int *n, int *k, double *A, int *lda, double *tau, double *work, int *lwork, int *info);
    void dorgqr_dlarfb0c2_(int *m, int *n, int *k, double *A, int *lda, double *tau, double *work, int *lwork, int *info);
    // LQ
    void dgelqf_(int *m, int *n, double *A, int *lda, double *tau, double *work, int *lwork, int *info);
    void dorglq_ref_(int *m, int *n, int *k, double *A, int *lda, double *tau, double *work, int *lwork, int *info);
    void dorglq_(int *m, int *n, int *k, double *A, int *lda, double *tau, double *work, int *lwork, int *info);
    void dorglq_dlarfb0c2_(int *m, int *n, int *k, double *A, int *lda, double *tau, double *work, int *lwork, int *info);
    // QL
    void dgeqlf_(int *m, int *n, double *A, int *lda, double *tau, double *work, int *lwork, int *info);
    void dorgql_ref_(int *m, int *n, int *k, double *A, int *lda, double *tau, double *work, int *lwork, int *info);
    void dorgql_(int *m, int *n, int *k, double *A, int *lda, double *tau, double *work, int *lwork, int *info);
    void dorgql_dlarfb0c2_(int *m, int *n, int *k, double *A, int *lda, double *tau, double *work, int *lwork, int *info);
    // RQ
    void dgerqf_(int *m, int *n, double *A, int *lda, double *tau, double *work, int *lwork, int *info);
    void dorgrq_ref_(int *m, int *n, int *k, double *A, int *lda, double *tau, double *work, int *lwork, int *info);
    void dorgrq_(int *m, int *n, int *k, double *A, int *lda, double *tau, double *work, int *lwork, int *info);
    void dorgrq_dlarfb0c2_(int *m, int *n, int *k, double *A, int *lda, double *tau, double *work, int *lwork, int *info);
#endif
