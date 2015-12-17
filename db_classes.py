import pdb
import MySQLdb
import math
import pickle
class DAO(object):
    db = None
    def __init__(self):
        #self.db = sqlite3.connect(DATABASE)
        self.db = MySQLdb.connect(host="localhost", # your host, usually localhost
                    user="root", # your username
                    passwd="genome", # your password
                    db="ES_OUTPUT") # name of the data base

    def exec_fetch_SQL(self,sql_template):
        cur = self.db.cursor()
        cur.execute(sql_template)
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
        sql_template=('select id from Geg where GWAS="' + GWAS + 
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
    def do_insert_gene_p_q(self,pval,qval,Geg_id):
    
        sql_template = 'insert into gene_p_q(pval,qval,Geg_id) values(%s,%s,%s)'
        cur = self.db.cursor()
        try:
            cur.execute(sql_template,(pval,qval,Geg_id))
        except sqlite3.ProgrammingError:
            pdb.set_trace()
            a = 1        
        self.db.commit()
        return cur.lastrowid
     
    def insert_gene_p_q(self,gene_p_q):
        Geg_id = self.fetch_or_insert_Geg(gene_p_q.GWAS,gene_p_q.eQTL,gene_p_q.gene)
        sql_template = 'select id,pval,qval from gene_p_q where Geg_id=' + str(Geg_id) + ';'
        cur = self.db.cursor()
        cur.execute(sql_template)
        rows = cur.fetchall()
        if len(rows) > 1:   # there is an unique constraint on (GWAS,eQTL,gene)
            print 'table gene_p_q is contaminated!'
            pdb.set_trace()
            exit(-1)
        elif len(rows) == 1 and rows[0][1] == gene_p_q.pval and rows[0][2] == gene_p_q.qval:
            return rows[0][0],Geg_id
        else:   # len(rows) == 0
            return self.do_insert_gene_p_q(gene_p_q.pval,gene_p_q.qval,Geg_id),Geg_id  
####################################
    def do_insert_single_GSNP_eSNP_Gpval(self,GSNP,eSNP,Gpval,Geg_id):
        sql_template= ('insert into GSNP_eSNP_Gpval (GSNP,eSNP,Gpval,Geg_id) values(%s,%s,%s,%s)')
        cur = self.db.cursor()
        try:
            cur.execute(sql_template,(GSNP,eSNP,Gpval,Geg_id))
        except sqlite3.ProgrammingError:
            pdb.set_trace()
            a = 1        
        self.db.commit()

    def insert_single_GSNP_eSNP_Gpval(self,GSNP,eSNP,Gpval,Geg_id):
        sql_template = (' select id, Gpval from GSNP_eSNP_Gpval where Geg_id=' + str(Geg_id) + 
                        ' and GSNP="' + GSNP + '"' +
                        ' and eSNP="' + eSNP + '"' +
                        ';')   
        cur = self.db.cursor() 
        cur.execute(sql_template)
        rows = cur.fetchall()
        if len(rows) >  1:   # there is an unique constraint on (Geg_id,SNP_Name)
            print 'GSNP_eSNP_Gpval is contaminated'
            exit(-1)
        elif len(rows) == 1 and rows[0][1] == Gpval:
            return rows[0][0]
        else: #len(rows) == 0
            return self.do_insert_single_GSNP_eSNP_Gpval(GSNP,eSNP,Gpval,Geg_id)

    def insert_gene_SNPs(self,gene_SNPs,Geg_id):
        #Geg_id = self.fetch_or_insert_Geg(gene_SNPs.GWAS,gene_SNPs.eQTL,gene_SNPs.gene)
        for t in gene_SNPs.GSNP_eSNP_Gpvals:
            GSNP = t[0]
            eSNP = t[1]
            Gpval= t[2]
            self.insert_single_GSNP_eSNP_Gpval(GSNP,eSNP,Gpval,Geg_id)


#pair manipulation   
    def fetch_gene_p_q_by_GWAS_eQTL(self,GWAS,eQTL):
        sql_template = ('select Geg.GWAS,Geg.eQTL,Geg.gene,gene_p_q.pval,gene_p_q.qval from Geg,gene_p_q'
        ' where Geg.id = gene_p_q.Geg_id'
        ' and Geg.GWAS="' + GWAS +'"'
        ' and Geg.eQTL="' + eQTL + '";' )
        cur = self.db.cursor()
        cur.execute(sql_template)
        rows = cur.fetchall()
        return list(rows) 

    def get_lowest_30_genes_for_all_pairs(self,result_lists):
        individual_genes = []
        for result in result_lists:
            result.sort(key=lambda x:float(x[3]))    
            rows_with_lowest_30_pval = result[:30]
            curr_gene = [x[2] for x in rows_with_lowest_30_pval] 
            individual_genes.append(set(curr_gene))
        comm_genes = set.union(*individual_genes)
        comm_gene_names = list(comm_genes)
        comm_gene_names.sort() 
        return comm_gene_names

    def filter_result_dict_by_lowest_30_genes_for_each_pair(self,result_dict):
        filtered_gene_names = self.get_lowest_30_genes_for_all_pairs(result_dict.values())
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
            term_relatives += [term,term.title(),term.lower(),term.upper()]
            term_relatives += [term + "'s", term + "'" ]
            term_relatives += [term + 's', term + 'es']
        return term_relatives

    def fetch_pair_gene(self,GWAS_list,eQTL_list):
        GWASs_disease_term = GWAS_list.strip().split()
        GWASs_disease_term = self.gen_term_relatives(GWASs_disease_term) 
        disease_GWAS_tuples_list = pickle.load(open(os.getcwd() + '/disease_GWAS_tuples_list.pickle','r'))
        GWASs = []
        # for each search_term, go over all disease_gwas tuple
        for disease in GWASs_disease_term:
            for disease_GWAS_tuple in disease_GWAS_tuples_list:
                if disease in disease_GWAS_tuple[0]:
                    GWASs.append(disease_GWAS_tuple[1]) 

        eQTLs = eQTL_list.strip().split()
        result_dict = {}
        for i in range(len(GWASs)):
            for j in range(len(eQTLs)):
                GWAS = GWASs[i]
                eQTL = eQTLs[j]
                result = self.fetch_gene_p_q_by_GWAS_eQTL(GWAS,eQTL)
                if len(result) > 0:
                    result_dict[GWAS + '---' + eQTL] = result
        filtered_dict,filtered_gene_names = self.filter_result_dict_by_lowest_30_genes_for_each_pair(result_dict)        
        return filtered_dict,filtered_gene_names

#pair manipulation   
#detail manipulation
    def fetch_SNP_list_by_GWAS_eQTL_gene(self,GWAS,eQTL,gene):
        sql_template = ('select GSNP,eSNP,Gpval from Geg, GSNP_eSNP_Gpval'
                        ' where Geg.GWAS = "' + GWAS + '"'
                        ' and Geg.eQTL = "' + eQTL + '"'
                        ' and Geg.gene = "' + gene + '"'
                        ' and GSNP_eSNP_Gpval.Geg_id = Geg.id;'  
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
        all_SNPs_list.sort()
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
                    curr_patched_list.append(('dummy','dummy','-1'))  
            patched_dict[pair_name] = curr_patched_list

        return patched_dict,all_SNPs_list


    def fetch_pair_SNP(self,GWAS_list,eQTL_list,gene):
        GWASs = GWAS_list.strip().split()
        eQTLs = eQTL_list.strip().split()

        result_dict = {}
        for i in range(len(GWASs)):
            for j in range(len(eQTLs)):
                GWAS = GWASs[i]
                eQTL = eQTLs[j]
                result = self.fetch_SNP_list_by_GWAS_eQTL_gene(GWAS,eQTL,gene)
                #result = self.fetch_gene_p_q_by_GWAS_eQTL(GWAS,eQTL)
                if len(result) > 0:
                    result_dict[GWAS + '---' + eQTL] = result
        patched_dict,all_SNPs_list = self.patch_result_dict_by_all_SNPs(result_dict)        
       
        return patched_dict,all_SNPs_list


#detail manipulation
    def fetch_detail(self,GWAS,eQTL,gene):
        sql_template = ('select GSNP,eSNP,Gpval from Geg, GSNP_eSNP_Gpval'
                        ' where Geg.GWAS = "' + GWAS + '"'
                        ' and Geg.eQTL = "' + eQTL + '"'
                        ' and Geg.gene = "' + gene + '"'
                        ' and GSNP_eSNP_Gpval.Geg_id = Geg.id;'  
                        )
        list_detail = self.exec_fetch_SQL(sql_template)
        return list_detail 
