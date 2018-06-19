#*************************
#   Created on 30.04.2018
#   author: Bita Khalili
#*************************

import numpy as np
import pandas as pd
from scipy import stats
#import sys
import matplotlib.pyplot as plt
from copy import deepcopy



def stocsy(M,ppm,dist,sig):
    print("\n----------Computing correlation matrix------------")
    M=correlation(M)
    print('Correlation matrix: '+\
          '{:d}'.format(M.shape[0])+'x'+'{:d}'.format(M.shape[1]))
        
    print("\n----------Computing off-diagonal distance------------")
    D=distanceMatrix(M.shape[0],ppm)
    output=pd.DataFrame(D)
    output.to_csv('distMatrix.csv',index=False,header=None)
    print('Correlation matrix: '+\
        '{:d}'.format(D.shape[0])+'x'+'{:d}'.format(D.shape[1]))
          
    print("\n----------Cunstructing significant correlation  matrix------------")
    M=M.values
    output=pd.DataFrame(M)
    output.to_csv('CorrMatrix.csv',index=False,header=None)

    P=maxCorrPairs(deepcopy(M),D,dist,ppm)
    P=remove_Pairs(P,M,ppm,sig,dist)

    print('Number of feature pairs with significant correlation : '+\
        '{:d}'.format(P.shape[0]))
        
    return(M,P)

def correlation(M):
    C=M.corr(method='pearson',min_periods=1)
    return(C)

def distanceMatrix(n,ppm):
    D = np.matmul(np.diag(ppm),np.ones((n,n)))   #a matrix with repeatition of ppms as columns
    D = abs(D-np.transpose(D))
    return (D)

def maxCorrPairs(CorrMat,DistMat,dist,ppm):
    P = np.zeros((CorrMat.shape[0] ,2),dtype=int)                  #list of pairs indices with max corr
    flag_feat =np.logical_not(DistMat<dist)
    CorrMat[np.invert(flag_feat)]=0
    #output=pd.DataFrame(CorrMat)
    #output.to_csv('OffdiagCorrMat.csv')
    max_array=np.max(CorrMat,1)
    #print(max_array)
    #P=np.array(np.argmax(CorrMat,axis=0),np.arange(CorrMat.shape[0]))
    P[:,1]=np.argmax(CorrMat,axis=0)
    P[:,0]=range(CorrMat.shape[0])
    #output=pd.DataFrame(P)
    #output.insert(0,'val',max_array)
    #output.to_csv('Pairs.csv',index=False,header=None)
    P=P[np.argsort(max_array)][:]
    P=P[::-1][:]
    #output=pd.DataFrame(P)
    #output.insert(0,'val',np.sort(max_array)[::-1])
    #output.to_csv('SortedPairs.csv',index=False,header=None)
    return (P)
    
def remove_Pairs(P,M,ppm,sig,dist):    #removes repeated and unsignificant pairs
    removing_similar=np.ones(P.shape[0])
    v=[]
    for i in range(P.shape[0]):
        if M[P[i,0]][P[i,1]]<sig:
            removing_similar[i]=0
        for j in range(i+1,P.shape[0]):
            if (P[i][0]==P[j][1] and P[i][1]==P[j][0]): #add an or statement here for checking if any two pairs are within the dist threshold if yes the one with min correlation should be also eliminated
                removing_similar[j]=0
            elif (ppm[P[i][0]]<ppm[P[j][0]]+dist and ppm[P[i][0]]>ppm[P[j][0]]-dist and ppm[P[i][1]]<ppm[P[j][1]]+dist and ppm[P[i][1]]>ppm[P[j][1]]-dist) or (ppm[P[i][0]]<ppm[P[j][1]]+dist and ppm[P[i][0]]>ppm[P[j][1]]-dist and ppm[P[i][1]]<ppm[P[j][0]]+dist and ppm[P[i][1]]>ppm[P[j][0]]-dist):
                #if M[P[i,0]][P[i,1]]>sig and M[P[j,0]][P[j,1]]>sig:
                #v.append(np.array([ppm[P[i][0]],ppm[P[i][1]],M[P[i,0]][P[i,1]],ppm[P[j][0]],ppm[P[j][1]],M[P[j,0]][P[j,1]]]))
                if M[P[i,0]][P[i,1]]>M[P[j,0]][P[j,1]]:
                    removing_similar[j]=0
                else:
                    removing_similar[i]=0

    #output=pd.DataFrame(v)
    #output.to_csv('pairswoNeighboringSimilars.csv',index=False,header=None)
    removing_similar=removing_similar>0
    P=P[removing_similar]
    #output=pd.DataFrame(P)
    #output.to_csv('pairswosimilars.csv',index=False,header=None)
 
    return (P)


