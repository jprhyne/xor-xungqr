import numpy as np
def normalizedPerf(timeVecs, n, kVec):
    # timeVecs is an array of arrays
    dim1,dim2 = timeVecs.shape
    retVec = np.zeros((dim1,dim2))
    numerVec = n*kVec**2
    for i in range(dim1):
        # for each vec in timeVecs
        for j in range(dim2):
            retVec[i][j] = numerVec[j] / (10**9 * timeVecs[i][j])
    return retVec

def larfbPerf(m,n,k,time):
    tmp = 4.0*m*n*k - 2.0*n*k**2
    return tmp / (time*10**9)

def orgqrPerf(timeVecs, mVec, nVec, kVec):
    (numVecs,_) = timeVecs.shape
    retVecs = np.zeros(timeVecs.shape)
    for i in range(numVecs):
        for j in range(len(mVec)):
            m = mVec[j]
            n = nVec[j]
            k = kVec[j]
            retVecs[i,j] = (4*m*n*k - 2*(m+n)*k**2 + 4.0/3.0*k**3 + 3*n*k - m*k - k**2 - 4.0/3.0*k) / (10**9 * timeVecs[i,j])
    return retVecs

def orgkrPerf(timeVecs, mVec, nVec):
    (numVecs,_) = timeVecs.shape
    retVecs = np.zeros(timeVecs.shape)
    for i in range(numVecs):
        for j in range(len(mVec)):
            retVecs[i,j] = 4*mVec[j] * nVec[j] * nVec[j] / (10**9 * timeVecs[i,j])
    return retVecs