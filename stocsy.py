#*************************
#   Created on 30.04.2018
#   author: Bita Khalili
#*************************

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from copy import deepcopy
import sys



def stocsy(M,ppm_labels,sig):
    loop_flag=True
    pseudospec_cr_labels=[]
    labels=deepcopy(ppm_labels)
    feature_matrix=deepcopy(M)   #original data Matrix (M)
    correlation_matrix=correlation(feature_matrix)
    index_vector=np.arange(M.shape[1]) #to contain the original indices
    feature_pair_index=[[index_vector[i]] for i in range(len(index_vector))] # list of lists to contain all the indices of features that are correlated above threshold
    iter=0
    while loop_flag:
        print("\n----------Computing correlation matrix------------")
        C=correlation(M)
        print('Correlation matrix: '+\
            '{:d}'.format(C.shape[0])+'x'+'{:d}'.format(C.shape[1]))
          
        print("\n----------Cunstructing significant correlation  matrix------------")
        C=C.values  #Correlation matrix
        #output=pd.DataFrame(M)
        #output.to_csv('CorrMatrix.csv',index=False,header=None)

        P=maxCorrPairs(deepcopy(C)) # indices of the maximum element in the new correlation matrix (at each step)
        print(P[0])    
        print(P[1])
        print(labels[P[0][0]])
        print(labels[P[1][0]])
        x=['f',str(labels[P[1][0]]),'f',str(labels[P[0][0]])]
        #print(x)
        labels[P[1][0]]=''.join(item.replace('\'','') for item in x)
        labels=np.delete(labels,P[0][0])
        print(labels[P[1][0]])
        pseudospec_cr_labels.append(ppm_labels[P[1][0]])
        feature_pair_index[index_vector[P[1][0]]]
        feature_pair_index[index_vector[P[1][0]]]=feature_pair_index[index_vector[P[1][0]]]+feature_pair_index[index_vector[P[0][0]]]     #the original index in the original feature-feature matrix
        print(feature_pair_index[P[1][0]])
        #print(sum([feature_matrix[:,x] for x in feature_pair_index[index_vector[P[1][0]]]]))
        print(len(feature_pair_index[index_vector[P[1][0]]]))
        M[:,index_vector[P[1][0]]]=sum([feature_matrix[:,items] for items in feature_pair_index[index_vector[P[1][0]]]])/len(feature_pair_index[index_vector[P[1][0]]])  #P holds the indices of max correlation value and feature_matrix holds the original feature vectors 
        print([feature_matrix[0,items] for items in feature_pair_index[index_vector[P[1][0]]]])
        #print(M[0,P[1][0]])
        M=np.delete(M,P[0][0],1)

        print('-------------')
        print([ppm_labels[items] for items in feature_pair_index[index_vector[P[1][0]]]])
        print('------------')
        index_vector=np.delete(index_vector,P[0][0])
        
        if iter==10:
            loop_flag=False
        iter=iter+1
        if C[P]<sig:
            loop_flag=False
    CorMat=correlation(feature_matrix)
    CorMat=CorMat.values
    output=pd.DataFrame(pseudospec_cr_labels)
    output.to_csv("pseudospec_cr_labels.tsv",index=False,header=None,sep='\t')
    output=pd.DataFrame(feature_pair_index)
    output.to_csv("feature_pair_index.tsv",index=False,header=None,sep='\t')
    print('Correlation matrix: '+\
            '{:d}'.format(CorMat.shape[0])+'x'+'{:d}'.format(CorMat.shape[1]))    
    return(CorMat,feature_pair_index)

def correlation(M):
    M=pd.DataFrame(M)
    C=M.corr(method='pearson',min_periods=1)
    return(C)


def maxCorrPairs(CorrMat):
    flag_feat=np.logical_not(np.diag(np.diag(CorrMat)))
    CorrMat[np.invert(flag_feat)]=0
    max_cr=np.max(np.nanmax(CorrMat,1))
    argmax_cr=np.where(np.tril(CorrMat)==(max_cr))
    # print(argmax_cr)
    # print(CorrMat[argmax_cr])
    # sys.exit()

    return (argmax_cr)
    



