# Abstract
The LAPACK routine dorg2r is used to attain the full Q matrix in the QR decomposition, and is needed even when 
the blocked version is used. Using a recursive scheme, we improve the performance of the dorg2r in the case of 
having number of reflectors equal to the number of columns of Q. In addition, our scheme produces a way to compute
The first $k$ columns of Q only, which can be useful when coupled with a more efficient scheme to produce the 
remaining columns needed. We also propose a new version of dorgqr where we also see performance increases for 
modestly sized inputs against both reference LAPACK and tested optimized BLAS routines.

# Notes for reference
* Mention dorg2r, dorgqr
* Mention the difference in outputs of our scheme vs dorg2r
* Discuss that we get performance gain
* Discuss the connection between our dorgqr work and dorg2r (Q_1 and Q_2 computations)
