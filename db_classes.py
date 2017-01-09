import pdb
import MySQLdb
import math
import pickle
import os
import time
import multiprocessing
import copy_reg,types
import psutil
import copy
from beans import Chrom_fields
class DAO(object):
    db = None
    db_hg19 = None
    db_GWAS = None
    display_name_GWAS_eQTL_tuple_dict = None

    EMPTY_CHROM_POS = 3200000000 
    EMPTY_CHR = 'CHROMOSOME NOT FOUND'
    DUMMY_SNP_NAME = 'dummy'
    def __init__(self,web_disease_list,web_eQTL_list):
        self.db = MySQLdb.connect(host="genomesvr2", # your host, usually localhost
                    user="es", # your username
                    passwd="detective", # your password
                    db="ES_OUTPUT") # name of the data base
 
        self.db_hg19 = MySQLdb.connect(host="genomesvr2", # your host, usually localhost
                    user="es",
                    passwd="detective",
                    db="hg19")

        #self.db = MySQLdb.connect(host="localhost", # your host, usually localhost
        #            user="root", # your username
        #            passwd="genome", # your password
        #            db="ES_OUTPUT") # name of the data base
        
        #self.db_hg19 = MySQLdb.connect(host="localhost",
        #            user="root",
        #            passwd="genome",
        #            db="hg19")

        #self.db_GWAS = MySQLdb.connect(host="localhost",
        #            user="root",
        #            passwd="genome",
        #            db="GWAS")
        if web_disease_list is not None and web_eQTL_list is not None:
            self.gen_display_name_GWAS_eQTL_tuple_dict(web_disease_list,web_eQTL_list)

    def gen_display_name_GWAS_eQTL_tuple_dict(self,web_disease_list,web_eQTL_list):
        Merged_name = 'Merged_08212015_pruned_LD02'

        GWASs,GWAS_disease_dict = self.gen_GWASs_from_web_disease_list(web_disease_list) 
        eQTLs = web_eQTL_list.strip().split()
        
        display_name_GWAS_eQTL_tuple_dict = {}
        if 'merged_pickle' in eQTLs:
            eQTLs.remove('merged_pickle')
            eQTLs.append(Merged_name)
        
        for i in range(len(GWASs)):
            for j in range(len(eQTLs)):
                GWAS = GWASs[i]
                eQTL = eQTLs[j]
                display_name = self.gen_display_name_from_GWAS_eQTL(GWAS_disease_dict,GWAS,eQTL,Merged_name)
                display_name_GWAS_eQTL_tuple_dict[display_name] = (GWAS,eQTL)

        self.display_name_GWAS_eQTL_tuple_dict = display_name_GWAS_eQTL_tuple_dict


    def exec_fetch_SQL(self,sql_template):
        cur = self.db.cursor()
        cur.execute(sql_template)
        rows = cur.fetchall()
        return list(rows) 

    def exec_fetch_SQL_hg19(self,sql_template):
        cur = self.db_hg19.cursor()
        cur.execute(sql_template)
        rows = cur.fetchall()
        return list(rows) 


    def exec_fetch_SQL_GWAS(self,sql_template):
        
        cur = self.db_GWAS.cursor()
        try:
            cur.execute(sql_template)
        except ProgrammingError as err:
            if err.errno == errorcode.ER_SYNTAX_ERROR:
                print("Check your syntax!")
                return None
            else:
                print("Error: {}".format(err))
                return None
        rows = cur.fetchall()
        return list(rows) 



    def do_insert_Geg(self,GWAS,eQTL,gene):
        sql_template = 'insert into Geg(GWAS,eQTL,gene) values(%s,%s,%s)'
        cur = self.db.cursor()
        try:
            cur.execute(sql_template,(GWAS,eQTL,gene))
        except sqlite3.ProgrammingError:
            pdb.set_trace()
            a = 1        
        self.db.commit()
        return cur.lastrowid

    def fetch_or_insert_Geg(self,GWAS,eQTL,gene):
        #sql_template=('select id from Geg where GWAS="' + GWAS + 
        #              '" and eQTL="' + eQTL + 
        #              '" and gene="'+ gene + 
        #              '";')
        sql_template=('select id from gene_fields where GWAS="' + GWAS + 
                      '" and eQTL="' + eQTL + 
                      '" and gene="'+ gene + 
                      '";')
        

        cur = self.db.cursor()
        cur.execute(sql_template)
        rows = cur.fetchall()
        if len(rows)> 1: 
            print 'table Geg is contaminated!' 
            exit(-1)
        elif len(rows) == 1:
            pdb.set_trace()
            return rows[0][0]   #the id of table Geg
        else:
            return self.do_insert_Geg(GWAS,eQTL,gene)           # this is a new GWAS_eQTL_gene combination

#####################################
    def do_insert_single_SNP_fields(self,GSNP,eSNP,Gpval,gene_fields_id):
        sql_template= ('insert into SNP_fields (GSNP,eSNP,Gpval,gene_fields_id) values(%s,%s,%s,%s)')
        cur = self.db.cursor()
        try:
            cur.execute(sql_template,(GSNP,eSNP,Gpval,gene_fields_id))
        except sqlite3.ProgrammingError:
            pdb.set_trace()
            a = 1        
        self.db.commit()

    def insert_single_SNP_fields(self,GSNP,eSNP,Gpval,gene_fields_id):
        sql_template = (' select id, Gpval from SNP_fields where gene_fields_id=' + str(gene_fields_id) + 
                        ' and GSNP="' + GSNP + '"' +
                        ' and eSNP="' + eSNP + '"' +
                        ';')   
        cur = self.db.cursor() 
        cur.execute(sql_template)
        rows = cur.fetchall()
        if len(rows) >  1:   # there is an unique constraint on (Geg_id,SNP_Name)
            print 'SNP_fields is contaminated'
            exit(-1)
        elif len(rows) == 1 and rows[0][1] == Gpval:
            return rows[0][0]
        else: #len(rows) == 0
            return self.do_insert_single_SNP_fields(GSNP,eSNP,Gpval,gene_fields_id)

    def insert_gene_SNPs(self,gene_SNPs,gene_fields_id):
        #Geg_id = self.fetch_or_insert_Geg(gene_SNPs.GWAS,gene_SNPs.eQTL,gene_SNPs.gene)
        for t in gene_SNPs.SNP_fields:
            GSNP = t[0]
            eSNP = t[1]
            Gpval= t[2]
            self.insert_single_SNP_fields(GSNP,eSNP,Gpval,gene_fields_id)


#pair manipulation 


    def fetch_gene_p_q_by_GWAS_Merged(self,GWAS):
        LD02 = self.fetch_gene_p_q_by_GWAS_eQTL(GWAS,'Merged_08212015_pruned_LD02')
        
        return LD02
        # uncomment this will choose smaller p-val out of LD85 and LD02
        if False: 
            LD85 = self.fetch_gene_p_q_by_GWAS_eQTL(GWAS,'Merged_08212015_pruned_LD85')
            #merge LD85 and LD02, keep the one with more significant p-value
            all_key = set()
            LD85_dict = {}
            LD02_dict = {}
            #populate dictionary for both 85 and 02
            for pair in LD85:
                gene_LD85 = pair[2]
                all_key.add(gene_LD85)
                LD85_dict[gene_LD85] = pair
            for pair in LD02:
                gene_LD02 = pair[2]
                all_key.add(gene_LD02)
                LD02_dict[gene_LD02] = pair
            #pick the one not NULL or has a more significant p-value
            Merged = [] 
            for key in all_key:
                if  key not in LD85_dict or float(LD85_dict[key][3]) < 0 : 
                    Merged.append(LD02_dict[key]) 
                    continue
                elif key not in LD02_dict or float(LD02_dict[key][3]) < 0 :
                    Merged.append(LD85_dict[key])
                    continue
                else:
                    pair_LD85 = LD85_dict[key]
                    pair_LD02 = LD02_dict[key]
                    pval_LD85 = float(pair_LD85[3])
                    pval_LD02 = float(pair_LD02[3])
                    if pval_LD85 < pval_LD02:
                        Merged.append(pair_LD85)
                    else:
                        Merged.append(pair_LD02)
            pdb.set_trace() 
            return Merged 

  
    def fetch_gene_p_q_by_GWAS_eQTL(self,GWAS,eQTL):
        start_time = time.time()
        #sql_template = ('select Geg.GWAS,Geg.eQTL,Geg.gene,gene_p_q.pval,gene_p_q.qval from Geg,gene_p_q'
        #' where Geg.id = gene_p_q.Geg_id'
        #' and Geg.GWAS="' + GWAS +'"'
        #' and Geg.eQTL="' + eQTL + '";' )
        sql_template = (' select GWAS,eQTL,gene,pval,qval from gene_fields'
                        ' where GWAS="' + GWAS + '"'
                        ' and eQTL="' + eQTL + '";')
        cur = self.db.cursor()
        cur.execute(sql_template)
        rows = cur.fetchall()
        #print('#### fetching gene_p_qs for {} and {} takes {}'.format(GWAS,eQTL,(time.time() - start_time)))
        return list(rows) 

    def pool_worker_fetch_gene_p_q_by_GWAS_eQTL(self,pair):
        #start_time = time.time()
        #print '$$$$$$$$$$$inside of pool_worker'
        GWAS = pair[0]
        eQTL = pair[1]
        #sql_template = ('select Geg.GWAS,Geg.eQTL,Geg.gene,gene_p_q.pval,gene_p_q.qval from Geg,gene_p_q'
        #' where Geg.id = gene_p_q.Geg_id'
        #' and Geg.GWAS="' + GWAS +'"'
        #' and Geg.eQTL="' + eQTL + '";' )
        
        sql_template = ('select GWAS,eQTL,gene,pval,qval from gene_fields'
        ' where GWAS="' + GWAS +'"'
        ' and   eQTL="' + eQTL + '";' )




        cur = self.db.cursor()
        cur.execute(sql_template)
        rows = cur.fetchall()
        #print('### fetching gene_p_qs for {} and {} takes {}'.format(GWAS,eQTL,(time.time() - start_time)))
        return list(rows) 
    # get the union of all gene_names for all GWAS-eQTL pairs
    def get_all_gene_names(self,result_dict):
        all_gene_names = set()
        for result in result_dict:
            curr_tuple_list = result_dict[result]
            curr_gene_names = [x[2] for x in curr_tuple_list]
            all_gene_names = all_gene_names.union(set(curr_gene_names))
        all_gene_names = list(all_gene_names)
        all_gene_names.sort()
        return all_gene_names

    def get_lowest_n_genes_for_all_pairs(self,result_lists,num_genes_per_pair):
        individual_genes = []
        for result in result_lists:
            result.sort(key=lambda x:float(x[3]))    
            rows_with_lowest_n_pval = result[:num_genes_per_pair]
            curr_gene = [x[2] for x in rows_with_lowest_n_pval] 
            individual_genes.append(set(curr_gene))
	    comm_genes = set.union(*individual_genes)
        comm_gene_names = list(comm_genes)
        comm_gene_names.sort() 
        return comm_gene_names

    # this function will also patch each pair to the same lenth with 'dummy genes'
    def filter_result_dict_by_lowest_n_genes_for_each_pair(self,result_dict,num_genes_per_pair,consider_all_genes_in_database):
        filtered_gene_names = []
        if not consider_all_genes_in_database:
            filtered_gene_names = self.get_lowest_n_genes_for_all_pairs(result_dict.values(),num_genes_per_pair)
        else:
            filtered_gene_names = self.get_all_gene_names(result_dict)
        filtered_dict = {}
        for pair_name in result_dict:
            result = result_dict[pair_name]
            curr_patched_list = []
            curr_result_dict = {}
            for gene_tuple in result:
                gene = gene_tuple[2]
                curr_result_dict[gene] = gene_tuple
            
            for gene in filtered_gene_names:
                if gene in curr_result_dict:
                    curr_patched_list.append(curr_result_dict[gene])
                else:
                    curr_patched_list.append(('dummy_GWAS','dummy_eQTL','dummy_gene','-1','-1'))  
            filtered_dict[pair_name] = curr_patched_list
        '''
        #lists contain the same set of genes, sorting each one by gene will render the lists in the same order in terms of genes
        for l in filtered_dict.values() :
            l.sort(key=lambda x:x[2])
        '''
        return filtered_dict,filtered_gene_names

    def gen_term_relatives(self,terms):
        term_relatives = []
        for term in terms:
            # term.title() will change chron's to Crohn'S, thus need to convert the 'S back to 's
            term_relatives += [term,term.title(),term.lower(),term.upper(), term.title().replace("'S","'s")]
            term_relatives += [term + "'s", term + "'" ]
            term_relatives += [term + 's', term + 'es']
        return term_relatives

    def gen_GWASs_from_web_disease_list(self,web_disease_list):
        web_diseases_term = web_disease_list.strip().split(',')
        web_diseases_term_found_dict = {}
        #disease_GWAS_dict = pickle.load(open(os.getcwd() + '/disease_GWAS_dict.pickle','r'))
        disease_GWAS_dict = self.build_disease_GWAS_dict() 
        GWASs = set([])
        # this dict is used to show which disease a certain GWAS belongs to
        GWAS_disease_dict = dict()
        # for each search_term, go over all disease_gwas tuple
        for disease in web_diseases_term:
            disease = disease.strip()
            if disease in disease_GWAS_dict:
                GWASs = GWASs.union(disease_GWAS_dict[disease])
                web_diseases_term_found_dict[disease] = True 
                matched_GWASs = disease_GWAS_dict[disease]
                for GWAS in matched_GWASs:
                    GWAS_disease_dict[GWAS] = disease
        
        # for any disease term not found a matching GWAS fwith exact match,
        # do fuzzy match
        for disease in web_diseases_term:
            # if not found in exact match
            if disease not in web_diseases_term_found_dict:
                disease_relatives = self.gen_term_relatives([disease])
                for relative in disease_relatives:
                    # for each relative, go over all keys
                    for full_disease_name in disease_GWAS_dict.keys():
                        if relative in full_disease_name:
                             GWASs = GWASs.union(disease_GWAS_dict[full_disease_name])
                             web_diseases_term_found_dict[disease] = True 
                             matched_GWASs = disease_GWAS_dict[full_disease_name]
                             for GWAS in matched_GWASs:
                                GWAS_disease_dict[GWAS] = relative

        return list(GWASs), GWAS_disease_dict

    def gen_display_name_from_GWAS_eQTL(self,GWAS_disease_dict,GWAS,eQTL,Merged_name):
        eQTL_tissue_dict = self.build_eQTL_tissue_dict()

        display_name = ''
        if eQTL == Merged_name:
            display_name = GWAS_disease_dict[GWAS] + '---Merged'  + '  (' + GWAS + '---Merged)'
        else:
            display_name = GWAS_disease_dict[GWAS] + '---' + eQTL_tissue_dict[eQTL] + '  (' + GWAS + '---' + eQTL + ')'
        return display_name 

    def fetch_pair_gene(self,web_disease_list,web_eQTL_list,web_num_genes_per_pair,consider_all_genes_in_database):
        Merged_name = 'Merged_08212015_pruned_LD02'

        GWASs,GWAS_disease_dict = self.gen_GWASs_from_web_disease_list(web_disease_list) 
        eQTLs = web_eQTL_list.strip().split()


        #contain_merged = False
        if 'merged_pickle' in eQTLs:
            #eQTLs.append('Merged_08212015_pruned_LD85')
            #contain_merged = True
            eQTLs.remove('merged_pickle')
            eQTLs.append(Merged_name)
        if (len(GWASs) == 0 or len(eQTLs) == 0): 
            return None,None,None
        num_genes_per_pair = 30
        try:
            num_genes_per_pair = int(web_num_genes_per_pair)
        except ValueError:
            print 'web_num_genes_per_pair is not an integer, use default value 30'


        #disease_GWAS_dict = pickle.load(open(os.getcwd() + '/disease_GWAS_dict.pickle','r'))
        #disease_GWAS_dict = self.build_disease_GWAS_dict()
        #eQTL_tissue_dict  = pickle.load(open(os.getcwd() + '/eQTL_tissue_dict.pickle','r'))
        #eQTL_tissue_dict = self.build_eQTL_tissue_dict()    

        start_time = time.time() 
        result_dict = {}
        for i in range(len(GWASs)):
            for j in range(len(eQTLs)):
                GWAS = GWASs[i]
                eQTL = eQTLs[j]
                result = self.fetch_gene_p_q_by_GWAS_eQTL(GWAS,eQTL)
                if len(result) > 0:
                    display_name = self.gen_display_name_from_GWAS_eQTL(GWAS_disease_dict,GWAS,eQTL,Merged_name)
                    result_dict[display_name] = result

        print("fetching all pairs using loop takes %s seconds" % (time.time() - start_time))
    
        if len(result_dict) == 0:
            return None,None,None
        filtered_dict,filtered_gene_names = self.filter_result_dict_by_lowest_n_genes_for_each_pair(result_dict,num_genes_per_pair,consider_all_genes_in_database)       
        #pickle.dump(filtered_dict,open('/genomesvr1/home/jgu1/WorkSpace/job_11_12/ES_web/longevity.pickle','wb+'))  
        #pickle.dump(filtered_gene_names,open('/genomesvr1/home/jgu1/WorkSpace/job_11_12/ES_web/gene_names.pickle','wb+'))  
        start_time = time.time()
        gene_descriptions = self.fetch_gene_descriptions_by_gene_names_in_memory(filtered_gene_names)
        #self.fetch_gene_descriptions_by_gene_names_in_memory(filtered_gene_names)
        print("get_gene_description takes {} seconds".format(time.time() - start_time))
        return filtered_dict,filtered_gene_names,gene_descriptions
 
    # fetch the entire table into memory and do matchup for genes in meory
    # tried fetching description one at a time, but the overhead is too significant from python to MySQL
    def fetch_gene_descriptions_by_gene_names_in_memory(self,gene_names):
        descriptions = []
        gene_description_dict = {}
        sql_template = 'select gd_app_sym, gd_app_name from HUGO;'
        cur = self.db_hg19.cursor()
        cur.execute(sql_template)
        rows = cur.fetchall()
        for row in rows:
            gene = row[0]
            description = row[1]
            gene_description_dict[gene] = description

        for i_gene in range(len(gene_names)):
            gene = gene_names[i_gene]
            if gene in gene_description_dict:
                descriptions.append(gene_description_dict[gene])
            else:
                descriptions.append('no description found in hg19 database')

        return descriptions

#pair manipulation   
#detail manipulation
    def fetch_SNP_list_by_GWAS_eQTL_gene(self,GWAS,eQTL,gene):
        #sql_template = ('select GSNP,eSNP,Gpval,epval from Geg, SNP_fields'
        #                ' where Geg.GWAS = "' + GWAS + '"'
        #                ' and Geg.eQTL = "' + eQTL + '"'
        #                ' and Geg.gene = "' + gene + '"'
        #                ' and SNP_fields.Geg_id = Geg.id;'  
        #                )
        sql_template = ('select GSNP,eSNP,Gpval,epval from gene_fields, SNP_fields'
                        ' where gene_fields.GWAS = "' + GWAS + '"'
                        ' and gene_fields.eQTL = "' + eQTL + '"'
                        ' and gene_fields.gene = "' + gene + '"'
                        ' and SNP_fields.gene_fields_id = gene_fields.id;'  
                        )

        list_detail = self.exec_fetch_SQL(sql_template)
        return list_detail        

    def get_comm_SNPs(self,result_lists):
        individual_SNPs = []
        for result in result_lists:
            curr_SNP = [x[0] for x in result]
            individual_SNPs.append(set(curr_SNP))
        comm_SNPs = set.intersection(*individual_SNPs)
        return comm_SNPs

    def filter_result_dict_by_comm_SNPs(self,result_dict):
        comm_SNPs = self.get_comm_SNPs(result_dict.values())
        filtered_dict = {}
        for pair_name in result_dict:
            result = result_dict[pair_name]
            curr_new_list = []
            for row in result:
                if row[0] in comm_SNPs:
                    curr_new_list.append(row)
            filtered_dict[pair_name] = curr_new_list

        for l in filtered_dict.values() :
            l.sort(key=lambda x:x[0])

        return filtered_dict

    def get_all_SNPs(self,result_lists):
        individual_SNPs = []
        for result in result_lists:
            curr_SNP = [x[0] for x in result]
            individual_SNPs.append(set(curr_SNP))
        all_SNPs = set.union(*individual_SNPs)
        all_SNPs_list = list(all_SNPs)
        return all_SNPs_list

    def patch_result_dict_by_all_SNPs(self,result_dict):
        all_SNPs_list = self.get_all_SNPs(result_dict.values())
        patched_dict = {}
        for pair_name in result_dict:
            result = result_dict[pair_name]
            curr_patched_list = []
            curr_result_dict = {}
            for SNP_tuple in result:
                SNP = SNP_tuple[0]
                curr_result_dict[SNP] = SNP_tuple
            for SNP in all_SNPs_list:
                if SNP in curr_result_dict:
                    curr_patched_list.append(curr_result_dict[SNP])
                else:
                    curr_patched_list.append(('dummy','dummy','-1','-1'))  
            patched_dict[pair_name] = curr_patched_list

        return patched_dict,all_SNPs_list


    def fetch_pair_SNP(self,web_disease_list,web_eQTL_list,gene):
        Merged_name = 'Merged_08212015_pruned_LD02'
        #GWASs = GWAS_list.strip().split()
        GWASs,GWAS_disease_dict = self.gen_GWASs_from_web_disease_list(web_disease_list) 
        eQTLs = web_eQTL_list.strip().split()

        if 'merged_pickle' in eQTLs:
            eQTLs.remove('merged_pickle')
            eQTLs.append(Merged_name)


        result_dict = {}
        for i in range(len(GWASs)):
            for j in range(len(eQTLs)):
                GWAS = GWASs[i]
                eQTL = eQTLs[j]
                result = self.fetch_SNP_list_by_GWAS_eQTL_gene(GWAS,eQTL,gene)
                #result = self.fetch_gene_p_q_by_GWAS_eQTL(GWAS,eQTL)
                if len(result) > 0:
                    display_name = self.gen_display_name_from_GWAS_eQTL(GWAS_disease_dict,GWAS,eQTL,Merged_name)
                    result_dict[display_name] = result
        
        patched_dict,all_SNPs_list = self.patch_result_dict_by_all_SNPs(result_dict)        
       
        return patched_dict,all_SNPs_list

    def fetch_pair_SNP_raw(self,web_disease_list,web_eQTL_list,gene):
        Merged_name = 'Merged_08212015_pruned_LD02'
        GWASs,GWAS_disease_dict = self.gen_GWASs_from_web_disease_list(web_disease_list) 
        eQTLs = web_eQTL_list.strip().split()

        if 'merged_pickle' in eQTLs:
            eQTLs.remove('merged_pickle')
            eQTLs.append(Merged_name)


        result_dict = {}
        for i in range(len(GWASs)):
            for j in range(len(eQTLs)):
                GWAS = GWASs[i]
                eQTL = eQTLs[j]
                result = self.fetch_SNP_list_by_GWAS_eQTL_gene(GWAS,eQTL,gene)
                #result = self.fetch_gene_p_q_by_GWAS_eQTL(GWAS,eQTL)
                if len(result) > 0:
                    display_name = self.gen_display_name_from_GWAS_eQTL(GWAS_disease_dict,GWAS,eQTL,Merged_name)
                    result_dict[display_name] = result
        
        patched_dict,all_SNPs_list = self.patch_result_dict_by_all_SNPs(result_dict)        
       
        return patched_dict,all_SNPs_list




#detail manipulation
    def fetch_detail(self,GWAS,eQTL,gene):
        #sql_template = ('select GSNP,eSNP,Gpval from Geg, SNP_fields'
        #                ' where Geg.GWAS = "' + GWAS + '"'
        #                ' and Geg.eQTL = "' + eQTL + '"'
        #                ' and Geg.gene = "' + gene + '"'
        #                ' and SNP_fields.Geg_id = Geg.id;'  
        #                )
        sql_template = ('select GSNP,eSNP,Gpval from gene_fields, SNP_fields'
                        ' where gene_fields.GWAS = "' + GWAS + '"'
                        ' and gene_fields.eQTL = "' + eQTL + '"'
                        ' and gene_fields.gene = "' + gene + '"'
                        ' and SNP_fields.gene_fields_id = gene_fields.id;'  
                        )


        list_detail = self.exec_fetch_SQL(sql_template)
        return list_detail

    def build_disease_GWAS_dict(self):
        sql_template = 'select disease,GWAS from disease_GWAS;'
        disease_GWAS_list = self.exec_fetch_SQL(sql_template)
        disease_GWAS_dict = {}
        disease_GWAS_list.sort(key = lambda x: x[0])
        prev_disease = disease_GWAS_list[0][0]
        prev_GWAS_list = []
        for i_tuple in range(len(disease_GWAS_list)):
            curr_tuple = disease_GWAS_list[i_tuple]
            curr_disease = curr_tuple[0]
            curr_GWAS    = curr_tuple[1]
            if curr_disease != prev_disease:
                disease_GWAS_dict[prev_disease] = prev_GWAS_list
                prev_GWAS_list = []
                prev_disease = curr_disease
            prev_GWAS_list.append(curr_GWAS)

        disease_GWAS_dict[curr_disease] = prev_GWAS_list # at the end of iteration, add last disease into dict

        return disease_GWAS_dict

    def build_eQTL_tissue_dict(self):
        sql_template = 'select eQTL,tissue from eQTL_tissue;'
        eQTL_tissue_list = self.exec_fetch_SQL(sql_template)
        eQTL_tissue_dict = {}
        for row in eQTL_tissue_list:
            eQTL   = row[0]
            tissue = row[1] 
            eQTL_tissue_dict[eQTL] = tissue
        return eQTL_tissue_dict


    def Manhattan_get_all_available_GWAS_in_db(self):
        sql_template = 'show tables;'
        table_names = self.exec_fetch_SQL_GWAS(sql_template)
        available_GWASs = []
        for row in table_names:
            available_GWASs.append(row[0])
        return available_GWASs

    #   given a SNP list, return a dict {SNP_name : (chrom,abs_location)}
    def Manhattan_build_snp_location_dict(self,SNP_set,chrom_abs_dict):
        SNP_list_str = '( ' + ','.join(SNP_set) + ')'
        sql_template = 'select name,chrom,chromStart from snp138 where name in' + SNP_list_str + ';'
        rows = self.exec_fetch_SQL_hg19(sql_template) 
        snp_location_dict = {} 
        for row in rows:
            name = row[0]
            chrom = row[1]
            chromStart = row[2]
            #if name in snp_location_dict:
            #    print 'dup SNP in snp138 for ' + name
            abs_location = chromStart
            if chrom in chrom_abs_dict:
                abs_location = abs_location + chrom_abs_dict[chrom]
            else:
                continue # one SNP can have multiple rows in snp138, skip the rows with unidentified chr to prevent overwriting previously correct dict value
                
            snp_location_dict[name] = (chrom,abs_location)
        return snp_location_dict

    # generate a set of all SNPs in pair_SNP_dict. 
    # all the GWAS and eQTL SNPs for a particular gene across all pairs
    def Manhattan_gen_SNP_set(self,pair_SNP_dict):
        all_SNP_set = set()
        for pair, curr_SNP_list in pair_SNP_dict.iteritems():
            # each SNP is in the format of (GSNP,eSNP,Gpval,epval)
            for SNP in curr_SNP_list:
                GSNP = '"' + SNP[0] + '"' 
                eSNP = '"' + SNP[1] + '"' 
                all_SNP_set.add(GSNP)
                all_SNP_set.add(eSNP) 
        return all_SNP_set

    # the GWAS_tuple and eQTL_tuple used to draw black dots in Manhattan plot
    # format defined on 12.20.2016:
    # (GSNP_name,GSNP_chr,GSNP_abs,GSNP_pval)
    # (eSNP_name,eSNP_chr,eSNP_abs,eSNP_pval,gene_name,GSNP_name)
    def Manhattan_build_GWAS_tuple_and_eQTL_tuple(self,SNP_tuple):
        GSNP_name = SNP_tuple[0]
        eSNP_name = SNP_tuple[1]
        GSNP_pval = SNP_tuple[2]
        eSNP_pval = SNP_tuple[3]
        
        GSNP_location = SNP_tuple[4]
        GSNP_chr = GSNP_location[0]
        GSNP_abs = GSNP_location[1]
            
        eSNP_location = SNP_tuple[5]
        eSNP_chr = eSNP_location[0]
        eSNP_abs = eSNP_location[1]

        gene_name = SNP_tuple[6]
        
        aligned = True
        tagged = True
 
        GWAS_tuple = (GSNP_name,GSNP_chr,GSNP_abs,GSNP_pval)
        eQTL_tuple = (eSNP_name,eSNP_chr,eSNP_abs,eSNP_pval,aligned,tagged,gene_name,GSNP_name)        
        return GWAS_tuple,eQTL_tuple

    #return 3 dict:
    # first dict: 
        #key: GWAS_eQTL pairname
        #value: list of tuple in the format (GSNP_name,eSNP_name,GSNP_pval,eSNP_pval,(GSNP_chr,GSNP_abs),(eSNP_chr,eSNP_abs),gene_name)
    #second dict:
        #key: GWAS_eQTL pairname
        #value: list of tuple in the format (GSNP_name,GSNP_chr,GSNP_abs,GSNP_pval)
    #third dict:
        #key: GWAS_eQTL pairname
        #value: list of tuple in the format (eSNP_name,eSNP_chr,eSNP_abs,eSNP_pval,gene_name,GSNP_name)
    def Manhattan_enhance_SNP_tuple_with_abs_location(self,pair_SNP_dict,SNP_location_dict,gene):
        pair_SNP_dict_with_location = {}
        GWAS_SNPlist_dict = {}
        eQTL_SNPlist_dict_for_curr_gene = {}
        for pair_name, SNP_tuple_list in pair_SNP_dict.iteritems():
            SNP_tuple_list_with_location = []
            GWAS_SNP_tuple_list = []
            eQTL_SNP_tuple_list = []
            for SNP_tuple in SNP_tuple_list:
                GSNP_name = SNP_tuple[0]
                eSNP_name = SNP_tuple[1]
                if GSNP_name == self.DUMMY_SNP_NAME or eSNP_name == self.DUMMY_SNP_NAME:
                    continue    # no need for Manhattan plot to consider DUMMY SNPs
                GSNP_location = (self.EMPTY_CHR,self.EMPTY_CHROM_POS)
                eSNP_location = (self.EMPTY_CHR,self.EMPTY_CHROM_POS)
                try:
                    GSNP_location = SNP_location_dict[GSNP_name] 
                    eSNP_location = SNP_location_dict[eSNP_name]
                except KeyError:    # the SNP can be 'dummy'
                    a = 1

                # generate a new tuple in the format (GSNP_name, eSNP_name, GSNP_pval, eSNP_pval, GSNP_loc, eSNP_loc, gene)  
                SNP_tuple_with_location = copy.deepcopy(SNP_tuple)
                SNP_tuple_with_location = SNP_tuple_with_location + (GSNP_location, eSNP_location,gene)    
                
                # add this new tuple in new list
                SNP_tuple_list_with_location.append(SNP_tuple_with_location)
                
                # add formatted_tuple into lists
                GWAS_tuple,eQTL_tuple = self.Manhattan_build_GWAS_tuple_and_eQTL_tuple(SNP_tuple_with_location)
                GWAS_SNP_tuple_list.append(GWAS_tuple)
                eQTL_SNP_tuple_list.append(eQTL_tuple)

            pair_SNP_dict_with_location[pair_name] = SNP_tuple_list_with_location
            GWAS_SNPlist_dict[pair_name] = GWAS_SNP_tuple_list
            eQTL_SNPlist_dict_for_curr_gene[pair_name] = eQTL_SNP_tuple_list
        return pair_SNP_dict_with_location,GWAS_SNPlist_dict,eQTL_SNPlist_dict_for_curr_gene

    def Manhattan_enhance_pair_SNP_dict_with_location(self, pair_SNP_dict,gene):
        chrom_abs_dict = self.Manhattan_gen_chrom_abs_dict()
        all_SNP_set = self.Manhattan_gen_SNP_set(pair_SNP_dict)        
        SNP_location_dict = self.Manhattan_build_snp_location_dict(all_SNP_set,chrom_abs_dict) 
        return self.Manhattan_enhance_SNP_tuple_with_abs_location(pair_SNP_dict,SNP_location_dict,gene) 

    def Manhattan_gen_gene_location_dict(self,genes):
        chrom_abs_dict = self.Manhattan_gen_chrom_abs_dict()
        gene_list = []
        for gene in genes:
            gene_list.append('"' + gene + '"')
        gene_list_str = '(' + ','.join(gene_list) + ')'
        sql_template = 'select gene,chrom,chromStart,chromEnd from gene_location where gene in ' + gene_list_str + ';'
        rows = self.exec_fetch_SQL(sql_template) 
        gene_location_dict_available_in_db = {} #some gene may not in db
        for row in rows:
            gene        = row[0]
            chrom       = row[1]
            chromStart  = row[2]
            chromEnd    = row[3]
            if chromStart is None:
                print 'gene "' + gene + '" has chromStart NULL'
                chromStart = 0
            if chromEnd is None:
                print 'gene "' + gene + '" has chromEnd NULL'
                chromEnd = 0 
            chrom_abs_start = chrom_abs_dict[chrom] + chromStart
            chrom_abs_end   = chrom_abs_dict[chrom] + chromEnd
            gene_location_dict_available_in_db[gene] = (chrom_abs_start,chrom_abs_end)

        gene_location_dict = {} #for those gene not in db, assign(None,None) for location
        genes.sort()
        for i_gene in range(len(genes)):
            curr_gene = genes[i_gene]
            if curr_gene in gene_location_dict_available_in_db:
                gene_location_dict[curr_gene] = gene_location_dict_available_in_db[curr_gene]
            else:
                gene_location_dict[curr_gene] = (None,None)

        return gene_location_dict

    def Manhattan_gen_chrom_abs_dict(self):
        # build chrom_name -> abs_location diction
        Chrom_len_dict = Chrom_fields.Chrom_len_dict
        Chrom_name_list = ['chr1' ,'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                          'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20',
                          'chr21','chr22','chrX','chrY']
        chrom_abs_dict = {}
        for chrom_name in Chrom_name_list:
            abs_location = 0
            for i_chrom in range(Chrom_name_list.index(chrom_name)):
                curr_chrom_name = Chrom_name_list[i_chrom]
                abs_location = abs_location + Chrom_len_dict[curr_chrom_name]
            chrom_abs_dict[chrom_name] = abs_location
        # build chrom_name -> abs_location diction
        return chrom_abs_dict


    def Manhattan_add_abs_location(self,pair_SNP_dict):
        chrom_abs_dict = self.Manhattan_gen_chrom_abs_dict()


    def Manhattan_gen_chrom_starts(self):
        chrom_abs_dict = self.Manhattan_gen_chrom_abs_dict()       
        chrom_starts = chrom_abs_dict.values()
        chrom_starts.sort()
        return chrom_starts    

    #Manhattan_SNP_fields in the format of  (GSNP_name,abs_location,GSNP_pval,gene,chrom) 
    def Manhattan_build_Manhattan_SNP_fields_list_dict(self,pair_SNP_dict,gene):
        start_time = time.time()
        Chrom_len_dict = Chrom_fields.Chrom_len_dict
        Chrom_name_list = ['chr1' ,'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                          'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20',
                          'chr21','chr22','chrX','chrY']
        #available_GWASs = self.Manhattan_get_all_available_GWAS_in_db()
        Manhattan_SNP_fields_list_dict = {}
        for pair in pair_SNP_dict:
            GWAS_eQTL = self.display_name_GWAS_eQTL_tuple_dict[pair]
            GWAS = GWAS_eQTL[0]
            #if GWAS not in available_GWASs:
                #pdb.set_trace()
                #a = 1
            #    continue
            SNP_list = pair_SNP_dict[pair]
            curr_SNPlist = []
            for SNP_tuple in SNP_list:
                GSNP_name = SNP_tuple[0]
                GSNP_pval = SNP_tuple[2]
                #abs_location,chrom = self.Manhattan_get_SNP_abs_location_chrom(GSNP_name,GWAS,Chrom_len_dict,Chrom_name_list)
                #Manhattan_SNP_fields = (GSNP_name,abs_location,GSNP_pval,gene,chrom) 
                Manhattan_SNP_fields = (GSNP_name,GSNP_pval,gene) 
                curr_SNPlist.append(Manhattan_SNP_fields) 
            Manhattan_SNP_fields_list_dict[pair] = curr_SNPlist

        #pdb.set_trace()
        #a = 1   
        print 'Manhattan fetching SNPs for gene "' + gene + '" takes ' + str(time.time() - start_time) + ' seconds' 
        return Manhattan_SNP_fields_list_dict, Manhattan_SNP_fields_list_dict.keys()

    def Manhattan_gen_eQTL_SNPlist(self,location_pval_chrom_SNPlist_dict,genes):
        gene_str_list = []
        chrom_abs_dict = self.Manhattan_gen_chrom_abs_dict()
        for gene in genes:
            gene_str_list.append('"' + gene + '"')
        gene_str = '(' + ','.join(gene_str_list) + ')'
        
        eQTL_SNPlist_dict = {}
        for pair in location_pval_chrom_SNPlist_dict:
            GWAS_eQTL = self.display_name_GWAS_eQTL_tuple_dict[pair]
            eQTL = GWAS_eQTL[1]
            sql_template = 'select SNP,chromStart,pval,gene,chrom from eQTLs where eQTL = "' + eQTL + '" and gene in ' + gene_str + ';'
            rows = self.exec_fetch_SQL(sql_template)
            curr_eQTL_SNPlist = []
            for row in rows:
                SNP         = row[0]
                chromStart  = row[1]
                pval        = row[2]
                gene        = row[3]
                chrom       = row[4]

                chrom_location_whole = chrom_abs_dict[chrom] + chromStart                
                
                eQTL_SNPlist_ele = (SNP,chrom_location_whole,str(pval),eQTL,chrom)
                curr_eQTL_SNPlist.append(eQTL_SNPlist_ele)

            a = 1
            eQTL_SNPlist_dict[pair] = curr_eQTL_SNPlist
               
        return eQTL_SNPlist_dict 
    def Manhattan_get_all_eQTL_names(self):
        sql_template = "select distinct eQTL from eQTLs;"
        rows = self.exec_fetch_SQL(sql_template)
        all_eQTL_names = [row[0] for row in rows]
        return all_eQTL_names
 
    def remove_sub_clusters_from_same_disease(self,sub_clusters):
        sub_clusters_new = []
        for sub_cluster in sub_clusters:
            disease_set = set()
            row_comb = sub_cluster.row_comb
            for row in row_comb:
                GWAS_eQTL   = self.display_name_GWAS_eQTL_tuple_dict[row]
                GWAS        = GWAS_eQTL[0]
                disease     = self.get_disease_by_GWAS(GWAS)
                disease_set.add(disease)
            if len(disease_set) > 1:
                sub_clusters_new.append(sub_cluster)
            else:
                print 'remove one subcluster all belong to disease: ' + str(disease_set)
        return sub_clusters_new
    
    def get_disease_by_GWAS(self,GWAS):
        sql_template = 'select disease from disease_GWAS where GWAS = "'  + GWAS + '";'
        rows = self.exec_fetch_SQL(sql_template)
        if len(rows) != 1:
            print 'cannot find disease for GWAS ' + GWAS
            return None
        else:
            return rows[0][0]      
