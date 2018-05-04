#*************************
#   Created on 30.04.2018
#   author: Bita Khalili
#*************************

import numpy as np

def correlation(M):
    C=np.corrcoef(M,rowvar=0)
    return(C)

def distanceMatrix(n):
    D = np.indices((n,n))
    D = abs(D[0]-D[1])*1/np.sqrt(2)
    return (D)

def sigCorrMat(M,D,dist,sig):
    flag_feat =np.logical_and(M>sig , D>dist)
    N = []
    for j in range(flag_feat.shape[1]):
            if np.any(flag_feat[:][j]):   #np.any(flag_feat[i]):
                N.append(M[:][j])
    return (N)

def stocsy(M,dist,sig):
    print("\n----------Computing correlation matrix------------")
    M=correlation(M)
    print('Correlation matrix: '+\
              '{:d}'.format(M.shape[0])+'x'+'{:d}'.format(M.shape[1]))

    print("\n----------Computing off-diagonal distance------------")
    D=distanceMatrix(M.shape[0])
    print('Correlation matrix: '+\
          '{:d}'.format(D.shape[0])+'x'+'{:d}'.format(D.shape[1]))

    print("\n----------Cunstructing significant correlation  matrix------------")
    sigCorrM=sigCorrMat(M,D,dist,sig)
    print('Correlation matrix: '+\
          '{:d}'.format(len(sigCorrM))+'x'+'{:d}'.format(sigCorrM[0].shape[0]))
    return(sigCorrM)
