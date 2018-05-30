#*************************
#   Created on 30.04.2018
#   author: Bita Khalili
#*************************

import numpy as np
import pandas as pd
from scipy import stats
import sys
import matplotlib.pyplot as plt
from copy import deepcopy


def stocsy(M,ppm,dist,sig):
    print("\n----------Computing correlation matrix------------")
    M=correlation(M)
    print('Correlation matrix: '+\
          '{:d}'.format(M.shape[0])+'x'+'{:d}'.format(M.shape[1]))
        
    print("\n----------Computing off-diagonal distance------------")
    D=distanceMatrix(M.shape[0],ppm)
    print('Correlation matrix: '+\
        '{:d}'.format(D.shape[0])+'x'+'{:d}'.format(D.shape[1]))
          
    print("\n----------Cunstructing significant correlation  matrix------------")
    M=M.values
    #output=pd.DataFrame(M)
    #output.to_csv('CorrMatrix.csv',index=False,header=None)
    #print(M[0])
    #print('in stocsy.stocsy before maxCorrPairs')
    P=maxCorrPairs(deepcopy(M),D,dist,sig,ppm)
    #print(M[0])
    #print('in stocsy.stocsy after maxCorrPairs')
    print('Number of feature pairs with significant correlation : '+\
        '{:d}'.format(P.shape[0]))
    pseudospec=pseudospectrum(M,P)   #a list of pseusospecs each correspond to the order of pairs in P
    return(pseudospec)

def correlation(M):
    C=M.corr(method='pearson',min_periods=1)
    return(C)

def distanceMatrix(n,ppm):
    D = np.matmul(np.diag(ppm),np.ones((n,n)))   #a matrix with repeatition of ppms as columns
    D = abs(D-np.transpose(D))
    return (D)

def maxCorrPairs(CorrMat,DistMat,dist,sig,ppm):
    P = []                  #list of pairs with max corr
    for i in range(CorrMat.shape[0]):
        flag_feat =np.logical_and(CorrMat[i]>sig , DistMat[i]>dist)
        if np.any(flag_feat):   #if any of the off diagonal corrolation elements was significant
            featCorr=CorrMat[i]
            featCorr[np.invert(flag_feat)]=0     #to find the argmax among significant off-diagonal elements put to zero all the insignificant and diagonal elements of the featCorr
            P.append(np.array([i,np.argmax(featCorr)]))

            #print(np.array([i,np.argmax(featCorr)]))
    removing_similar=np.ones(len(P))
    for i in range(len(P)):
        for j in range(i+1,len(P)):
            if P[i][0]==P[j][1] and P[i][1]==P[j][0]:
                removing_similar[j]=0

    removing_similar=removing_similar>0
    #print(removing_similar)
    maxCorrPairMat=np.array(P)
    #print(maxCorrPairMat[:])
    #print(len(maxCorrPairMat))
    maxCorrPairMat=maxCorrPairMat[removing_similar]
    #print(maxCorrPairMat[:])
    #print(len(maxCorrPairMat))
    return (maxCorrPairMat)

def pseudospectrum(M,P):
    pseudospec=[]
    #print('*****in pseudospectrum*****')
    for i in range(len(P)):
        v=0.5*(M[P[i][0]]+M[P[i][1]])
        #print(v)
        pseudospec.append(v)
    return(pseudospec)
