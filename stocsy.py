#*************************
#   Created on 30.04.2018
#   author: Bita Khalili
#*************************

import numpy as np

def correlation(M):
    T=np.transpose(M)
    C=np.corrcoef(T)
    return(C)

def distanceMatrix(n,m):
    D = np.empty((n,m))
    for i in range(n):
        for j in range(m):
            D[i][j] = abs((i-0.5)/n-(j-0.5)/m)/np.sqrt(1/(n*n)+1/(m*m))
    return (D)

def sigCorrMat(M,D,dist,sig):
    flag_feat = np.zeros(M.shape[0],dtype=int)
    for j in range(D.shape[1]):
        for i in range(D.shape[0]):
            if(D[i][j] > dist):
                if(M[i][j] > sig):
                    flag_feat[j] = 1;
    N=[]
    for i in range(M.shape[0]):
        if flag_feat[i]:
            N.append(M[:][i])
    return (N)
