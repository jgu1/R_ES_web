import sys
import os
import pickle
import pdb
def Main():
    disease_GWAS_tuples_list = []
    disease_GWAS_dict={}
    count = 0
    for root, dirs, files in os.walk(GWAS_dir):
        for f in files:
            if f.endswith('.meta'):
                count = count  + 1    
                with open(root + '/' + f,'r' ) as metaFile:
                    for i,line in enumerate(metaFile):
                        if i == 1:
                            curr_disease = line.strip()                        
                            GWAS_name = f[:-5]
                            existing_GWAS_list = [] 
                            if curr_disease in disease_GWAS_dict:
                                existing_GWAS_list = disease_GWAS_dict[curr_disease]
                           
                            existing_GWAS_list.append(GWAS_name)
                            disease_GWAS_dict[curr_disease] = existing_GWAS_list
                                
                            #disease_GWAS_tuples_list.append((curr_disease,f[:-5]))
                            #disease_GWAS_dict[curr_disease] = f[:-5] #strip '.meta'
                            break
    pdb.set_trace()
    #return disease_GWAS_tuples_list
    return disease_GWAS_dict
if __name__=='__main__':
    GWAS_dir = sys.argv[1]
    GWAS_dir = '/home/jiashun/empirical_sherlock/data/GWAS'
    disease_GWAS_dict = Main()
    pik_file = open(os.getcwd() + '/disease_GWAS_dict.pickle','w')
    pickle.dump(disease_GWAS_dict,pik_file)
