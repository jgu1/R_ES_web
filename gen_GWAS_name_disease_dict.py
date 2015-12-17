import sys
import os
import pickle
import pdb
def Main():
    disease_GWAS_tuples_list = []
    for root, dirs, files in os.walk(GWAS_dir):
        for f in files:
            if f.endswith('.meta'):
                with open(root + '/' + f,'r' ) as metaFile:
                    for i,line in enumerate(metaFile):
                        if i == 1:
                            curr_disease = line.strip()                        
                            disease_GWAS_tuples_list.append((curr_disease,f[:-5]))
                            #disease_GWAS_dict[curr_disease] = f[:-5] #strip '.meta'
                            break
    pdb.set_trace()
    return disease_GWAS_tuples_list

if __name__=='__main__':
    GWAS_dir = sys.argv[1]
    disease_GWAS_tuples_list = Main()
    pik_file = open(os.getcwd() + '/disease_GWAS_tuples_list.pickle','w')
    pickle.dump(disease_GWAS_tuples_list,pik_file)
