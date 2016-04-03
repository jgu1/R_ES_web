import sys
import os
import pickle
import pdb
pik_file = 'eQTL_tissue_dict.pickle'
eQTL_dir = '/home/jiashun/empirical_sherlock/data/eQTL'
def Main():
    eQTL_tissue_dict={}
    for root, dirs, files in os.walk(eQTL_dir):
        for f in files:
            if f.endswith('.meta'):
                with open(root + '/' + f,'r' ) as metaFile:
                    for i,line in enumerate(metaFile):
                        if i == 1:
                            curr_tissue = line.strip()                        
                            eQTL_name = f[:-5]
                            eQTL_tissue_dict[eQTL_name] = curr_tissue
    pdb.set_trace()
    return eQTL_tissue_dict
if __name__=='__main__':
    eQTL_tissue_dict = Main()
    pik_file = open(os.getcwd() + '/' + pik_file,'w')
    pickle.dump(eQTL_tissue_dict,pik_file)
