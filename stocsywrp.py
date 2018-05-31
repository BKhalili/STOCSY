#*************************
#   Created on 30.04.2018
#   author: Bita Khalili
#*************************

import stocsy
import pandas as pd
from optparse import OptionParser
import sys
import math

def pseudospectrum(C,P,n):
    pseudospec=[]
    #print('*****in pseudospectrum*****')
    for i in range(len(P)):
        v=0.5*(C[P[i][0]]+C[P[i][1]])
        for j in range(v.shape[0]):
            v[j]=math.atanh(v[j])*math.sqrt(n-3)
        pseudospec.append(v)
    return(pseudospec)


def main():
    print("----------Asigning parameters------------")
    usage = "usage: %prog [optional] arg"
    parser = OptionParser(usage)
    parser.add_option('-i','--inputfile',dest='infile', type='string',default='data.csv')
    parser.add_option('-o','--outputfile',dest='outfile', type='string',default='stocsy.csv')
    parser.add_option('-d', '--dist',help='off-diagonal distance',dest='dist',type='float',default=0.1)
    parser.add_option('-s', '--significance',help='correlation significance',dest='sig',type='float',default=0.8)
    (options, args) = parser.parse_args()
    print("input file=",options.infile,"\noutput file=",options.outfile,"\noff-diagonal distance=",options.dist,"\ncorrolation significance=",options.sig)

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

    (C,P) = stocsy.stocsy(M,ppm,options.dist,options.sig)

    pseudospec = pseudospectrum(C,P,num_samples)   #a list of pseusospecs each correspond to the respective pairs in P

    print("\n----------writing output file------------")
    output=pd.DataFrame(pseudospec)
    output=output.transpose()
    output.insert(0,'shift',ppm)
    headers=[]
    for i in range(output.shape[1]):
        if(i==0):
            headers.append('shift')
        else:
            headers.append('z/C'+'{:03}'.format(i))
    
    output.to_csv(options.outfile,index=False,header=headers)

if __name__ == '__main__':
    main()

