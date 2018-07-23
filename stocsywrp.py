#*************************
#   Created on 30.04.2018
#   author: Bita Khalili
#*************************

import stocsy
import pandas as pd
from optparse import OptionParser
from copy import deepcopy
import numpy as np
import sys
import math
import os

def pseudospectrum(C,feature_pair_index,ppm_labels):    #M is the initial matrix with feature vectore and C is the complete correlation matrix
    ps_cr=[sum(C[items])/len(items) for items in feature_pair_index if len(items)>1]
    print(len(ps_cr))

    headers_cr=[]
    for items in feature_pair_index:
        if len(items)>1:
            print([ppm_labels[items[i]][0] for i in range(len(items)) if len(items)>1])
            headers_cr.append([''.join(ppm_labels[items[i]][0].replace('\'','')) for i in range(len(items))])  #''.join(ppm_labels[items].replace('\'','')) for items in feature_pair_index if len(items)>1
    print(len(headers_cr))
    # for i in range(len(P)):
    #     cr_v=0.5*(C[P[i][0]]+C[P[i][1]])
    #     pseudospec_cr.append(cr_v)
    #     featPseudospec[:,i]=0.5*(M[:,P[i][0]]+M[:,P[i][1]])
    #     z_v=np.empty_like(cr_v)
    #     for j in range(cr_v.shape[0]):
    #         z_v[j]=math.atanh(cr_v[j])*math.sqrt(n-3)
    #     pseudospec_z.append(z_v)
    return(ps_cr,headers_cr)

def rescale(pseudospec,r):
    pseudospec=[x/r for x in pseudospec]
    return(pseudospec)

def main():
    print("----------Asigning parameters------------")
    usage = "usage: %prog [optional] arg"
    parser = OptionParser(usage)
    parser.add_option('-i','--inputfolder',dest='infolder', type='string',default='/data/')
    parser.add_option('-o','--omit',dest='OmitNeigboringFlag',default=True)
    parser.add_option('-d', '--dist',help='off-diagonal distance',dest='dist',type='float',default=0.1)
    parser.add_option('-s', '--significance',help='correlation significance',dest='sig',type='float',default=0.48)
    parser.add_option('-x', '--numofps',help='number of pseudospectrums',dest='numofps',type='int',default=179)
    parser.add_option('-r', '--rescale',help='rescaling factor',dest='resc',type='float',default=1.0)
    (options, args) = parser.parse_args()
    print("input file=",options.infolder,"\nflag for omitting the feature-feature pairs belonging to same neigborhood=",options.OmitNeigboringFlag,"\noff-diagonal distance=",options.dist,"\ncorrolation significance=",options.sig,"\nz-score rescaling factor=",options.resc)
    source_dir=os.getcwd()
    if not os.path.exists(source_dir+options.infolder):
        print('Couldn not find the input folder!')
        sys.exit()
    
    for items in os.listdir(source_dir+options.infolder):
        if '.csv' in items:
            outdir='./data.out/ps.stocsy_'+items.split('.csv')[0]+'/'
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            print(source_dir+options.infolder+items)
            print("\n----------Reading input file------------")
            D = pd.read_csv(source_dir+options.infolder+items,header=None)
            M = D.loc[1:,1:]
            M = M.values
            num_samples=M.shape[0]
            ppm = D.loc[0,:]
            ppm = ppm.drop(labels=0)
            #ppm_labels=ppm.values.astype(str)
            ppm_labels=[]
            ppm=ppm.values.tolist()
            for i in range(len(ppm)):
            	ppm_labels.append(['f'+str(ppm[i])])
            print('data matrix loaded: '+\
                      '{:d}'.format(M.shape[0])+'x'+'{:d}'.format(M.shape[1]))
            (CorrMat,feature_pair_index)=stocsy.stocsy(deepcopy(M),deepcopy(ppm_labels),options.sig) #receives the correlation matrix of original feature data matrix and all the ppm labels that that might have been replaced by averaged in in stocsy
            #print('Number of feature pairs with significant correlation : '+\
                #'{:d}'.format(P.shape[0]))
            (pseudospec_cr,headers_cr)=pseudospectrum(CorrMat,feature_pair_index,deepcopy(ppm_labels))
            #pseudospec_z = pseudospectrum(C,num_samples)   #a list of pseusospecs each correspond to the respective pairs in P

            #pseudospec_z = rescale(pseudospec_z,options.resc)
            #change what comes below into a separated function write_output(pseudospec_z,ppm,pseudospec_cr,P)

            print("\n----------writing output files------------")
            #-------------Creating two files of pseoudospectrums from averaged correlation vectors and Fisher transformed of leading features
            #output_z=pd.DataFrame(pseudospec_z)
            #output_z=output_z.transpose()
            #output_z.insert(0,'shift',ppm)
            output_cr=pd.DataFrame(pseudospec_cr)
            output_cr=output_cr.transpose()

            #fileName_z=outdir+'stocsy_z_rescaledCoeff_'+'{:03}'.format(options.resc)+'_CorrSig_'+'{:03}'.format(options.sig)+'.pseudospectrum.tsv'
            #output_z.to_csv(fileName_z,index=False,header=headers_z,sep='\t')
            fileName_cr=outdir+'stocsy_cr_CorrSig_'+'{:03}'.format(options.sig)+'.pseudospectrum.tsv'
            output_cr.to_csv(fileName_cr,index=False,header=headers_cr,sep='\t')

            #-------------Creating a parameters in and descreption files
            # DescFile=open(outdir+'description.tsv','w')
            # del headers_z[0]
            # tags=[x.replace('z/','') for x in headers_z]
            # for i in range(len(tags)):
            #     DescFile.write(tags[i]+'\t')
            #     for j in range(P.shape[1]):
            #         if ((ppm[P[i,j]]*1e4%100)<50):
            #             DescFile.write('f'+str(round(ppm[P[i,j]],2))+'p')
            #         elif ((ppm[P[i,j]]*1e4%100)>50):
            #             DescFile.write('f'+str(round(ppm[P[i,j]],2))+'m')
            #         if j==0:
            #             DescFile.write('-')
            #     DescFile.write('\n')
            # DescFile.close()
            # paramFile=open(outdir+'parameters.in.tsv','w')
            # paramFile.write('reference'+'\tumdb\n')
            # paramFile.write('dsh'+'\t0.025\n')
            # paramFile.write('significant'+'\t5\n')
            # paramFile.write('suggestive'+'\t3\n')
            # paramFile.write('nshow'+'\t16\n')
            # paramFile.write('variant'+'\tpm\n')
            # paramFile.write('n_permutation'+'\t999\n')
            # paramFile.write('samplesize'+'\t'+str(num_samples)+'\n')
            # paramFile.close()
            # #-------------Creating a file of pseoudospectrums from averaged leading feature vectors
            # output=pd.DataFrame(featPseudospec)
            # headers_featPseudospec=[x.replace('z','f') for x in headers_z]
            # fileName=outdir+'featPseudospec'+'_CorrSig_'+'{:03}'.format(options.sig)+'.csv'
            # output.to_csv(fileName,index=False,header=headers_featPseudospec)
            # #---------------producing a table with f-f labels and their leading corresponding correlation
            # CorrV=np.ones((1,len(P)))
            # for i in range(len(P)):
            #     CorrV[0,i]=C[P[i][0],P[i][1]]

            # output=pd.DataFrame(CorrV)
            # headers_corr=[x.replace('z/','Corr/') for x in headers_z]
            # output.to_csv(outdir+'TableOfSortedHighlyCorrFeats.csv',index=False,header=headers_corr,sep='\t')



if __name__ == '__main__':
    main()

