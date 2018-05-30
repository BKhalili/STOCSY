#*************************
#   Created on 30.04.2018
#   author: Bita Khalili
#*************************

import stocsy
import pandas as pd
from optparse import OptionParser
import sys

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
    ppm = D.loc[0,:]
    ppm=ppm.drop(labels=0)
    ppm=ppm.values
    #fillna(0) filling any missing element by 0
    #M = M.values
    print('data matrix loaded: '+\
              '{:d}'.format(M.shape[0])+'x'+'{:d}'.format(M.shape[1]))

    pseudospectrums=stocsy.stocsy(M,ppm,options.dist,options.sig)

    print("\n----------writing output file------------")
    output=pd.DataFrame(pseudospectrums)
    output.to_csv(options.outfile,index=False,header=None)

if __name__ == '__main__':
    main()
