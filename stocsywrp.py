#*************************
#   Created on 30.04.2018
#   author: Bita Khalili
#*************************

import stocsy
import pandas as pd
from optparse import OptionParser
import numpy as np
import sys
import math

def pseudospectrum(M,C,P,n):    #M is the initial matrix with feature vectore and C is the complete correlation matrix
    pseudospec=[]
    featPseudospec=np.zeros((M.shape[0],len(P)))
    for i in range(len(P)):
        v=0.5*(C[P[i][0]]+C[P[i][1]])
        featPseudospec[:,i]=0.5*(M[:,P[i][0]]+M[:,P[i][1]])
        for j in range(v.shape[0]):
            v[j]=math.atanh(v[j])*math.sqrt(n-3)
        pseudospec.append(v)
    return(pseudospec,featPseudospec)

def rescale(pseudospec,r):
    pseudospec=[x/r for x in pseudospec]
    return(pseudospec)

def main():
    print("----------Asigning parameters------------")
    usage = "usage: %prog [optional] arg"
    parser = OptionParser(usage)
    parser.add_option('-i','--inputfile',dest='infile', type='string',default='./data/data.csv')
    parser.add_option('-o','--outputfile',dest='outfile', type='string',default='./data.out/')
    parser.add_option('-d', '--dist',help='off-diagonal distance',dest='dist',type='float',default=0.1)
    parser.add_option('-s', '--significance',help='correlation significance',dest='sig',type='float',default=0.8)
    parser.add_option('-r', '--rescale',help='rescaling factor',dest='r',type='float',default=1.0)
    (options, args) = parser.parse_args()
    print("input file=",options.infile,"output file=",options.outfile,"\noff-diagonal distance=",options.dist,"\ncorrolation significance=",options.sig,"\nz-score rescaling factor=",options.r)

    print("\n----------Reading input file------------")
    D = pd.read_csv(options.infile,header=None)
    M = D.loc[1:,1:]
    num_samples=M.shape[0]
    ppm = D.loc[0,:]
    ppm = ppm.drop(labels=0)
    ppm = ppm.values
    #fillna(0) filling any missing element by 0
    #M = M.values
    print('data matrix loaded: '+\
              '{:d}'.format(M.shape[0])+'x'+'{:d}'.format(M.shape[1]))

    (C,P) = stocsy.stocsy(M,ppm,options.dist,options.sig) #Complete correlation matrix and the sorted indices of leading correlated features
    M = M.values
    (pseudospec,featPseudospec) = pseudospectrum(M,C,P,num_samples)   #a list of pseusospecs each correspond to the respective pairs in P

    pseudospec = rescale(pseudospec,options.r)

    print("\n----------writing output file------------")
    #-------------Creating a file of pseoudospectrums from averaged correlation vectors of leading features
    output=pd.DataFrame(pseudospec)
    output=output.transpose()
    output.insert(0,'shift',ppm)
    headers=[]
    for i in range(output.shape[1]):
        if(i==0):
            headers.append('shift')
        else:
            headers.append('z/f'+'{:05}'.format(ppm[P[i-1,0]])+'f'+'{:05}'.format(ppm[P[i-1,1]]))

    fileName=options.outfile+'stocsy_rescaledCoeff_'+'{:03}'.format(options.r)+'_CorrSig_'+'{:03}'.format(options.sig)+'.csv'
    output.to_csv(fileName,index=False,header=headers)
    #-------------Creating a parameters in and descreption files
    DescFile=open(options.outfile+'description.tsv','w')
    del headers[0]
    tags=[x.replace('z/','') for x in headers]
    for i in range(len(tags)):
        DescFile.write(tags[i]+'\t')
        #+str(ppm[P[i,0]])+'-'+str(ppm[P[i,1]])+'\n')
        for j in range(P.shape[1]):
            if ((ppm[P[i,j]]*1e4%100)<50):
                DescFile.write('f'+str(round(ppm[P[i,j]],2))+'p')
            elif ((ppm[P[i,j]]*1e4%100)>50):
                DescFile.write('f'+str(round(ppm[P[i,j]],2))+'m')
            if j==0:
                DescFile.write('-')
        DescFile.write('\n')
    DescFile.close()
    paramFile=open(options.outfile+'parameters.in.tsv','w')
    paramFile.write('refrences'+'\tumdb\n')
    paramFile.write('dsh'+'\t0.03')
    paramFile.close()
    #-------------Creating a file of pseoudospectrums from averaged leading feature vectors
    output=pd.DataFrame(featPseudospec)
    headers_featPseudospec=[x.replace('z','f') for x in headers]
    fileName=options.outfile+'featPseudospec'+'_CorrSig_'+'{:03}'.format(options.sig)+'.csv'
    output.to_csv(fileName,index=False,header=headers_featPseudospec,sep='\t')
    #---------------producing a table with f-f labels and their leading corresponding correlation
    CorrV=np.ones((1,len(P)))
    for i in range(len(P)):
        CorrV[0,i]=C[P[i][0],P[i][1]]

    output=pd.DataFrame(CorrV)
    headers_corr=[x.replace('z/','Corr/') for x in headers]
    output.to_csv(options.outfile+'TableOfSortedHighlyCorrFeats.csv',index=False,header=headers_corr,sep='\t')



if __name__ == '__main__':
    main()

