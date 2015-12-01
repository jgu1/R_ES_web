import pdb
import MySQLdb
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
   
    def fetch_gene_p_q_by_GWAS_eQTL(self,GWAS,eQTL):
        sql_template = ('select Geg.GWAS,Geg.eQTL,Geg.gene,gene_p_q.pval,gene_p_q.qval from Geg,gene_p_q'
        ' where Geg.id = gene_p_q.Geg_id'
        ' and Geg.GWAS="' + GWAS +'"'
        ' and Geg.eQTL="' + eQTL + '";' )
        cur = self.db.cursor()
        cur.execute(sql_template)
        rows = cur.fetchall()
        return list(rows) 

    def get_comm_genes(self,result_lists):
        individual_genes = []
        for result in result_lists:
            curr_gene = [x[2] for x in result]
            individual_genes.append(set(curr_gene))
        comm_genes = set.intersection(*individual_genes)
        return comm_genes

    def filter_result_dict_by_comm_genes(self,result_dict):
        comm_genes = self.get_comm_genes(result_dict.values())
        filtered_dict = {}
        for pair_name in result_dict:
            result = result_dict[pair_name]
            curr_new_list = []
            for row in result:
                if row[2] in comm_genes:
                    curr_new_list.append(row)
            filtered_dict[pair_name] = curr_new_list
        #lists contain the same set of genes, sorting each one by gene will render the lists in the same order in terms of genes
        for l in filtered_dict.values() :
            l.sort(key=lambda x:x[2])

        return filtered_dict
    '''
    def fetch_all_gene_p_q(self):
        GWAS1 = 'Barrett_08'
        GWAS2 = 'Longevity_2014_Age85'
        GWAS3 = 'UCSF_AML'
        eQTL1 = 'merged_pickle'
        eQTL2 = 'merged_pickle'
        eQTL3 = 'merged_pickle'
        Barrett = self.fetch_gene_p_q_by_GWAS_eQTL(GWAS1,eQTL1)
        Longevity = self.fetch_gene_p_q_by_GWAS_eQTL(GWAS2,eQTL2)    
        UCSF_AML = self.fetch_gene_p_q_by_GWAS_eQTL(GWAS3,eQTL3)    
        result_lists = [Barrett,Longevity,UCSF_AML]
        new_lists = self.filter_result__by_comm_genes(result_lists)        

        gene_p_qs = {}
        gene_p_qs[GWAS1+eQTL1] = new_lists[0]
        gene_p_qs[GWAS2+eQTL2] = new_lists[1]
        gene_p_qs[GWAS3+eQTL3] = new_lists[2]
        return gene_p_qs
    '''
    def fetch_pair_gene(self,GWAS_list,eQTL_list):
        GWASs = GWAS_list.strip().split()
        eQTLs = eQTL_list.strip().split()

        result_dict = {}
        for i in range(len(GWASs)):
            for j in range(len(eQTLs)):
                GWAS = GWASs[i]
                eQTL = eQTLs[j]
                result = self.fetch_gene_p_q_by_GWAS_eQTL(GWAS,eQTL)
                if len(result) > 0:
                    result_dict[GWAS + eQTL] = result
        filtered_dict = self.filter_result_dict_by_comm_genes(result_dict)        
        return filtered_dict


    def fetch_detail(self,GWAS,eQTL,gene):
        sql_template = ('select GSNP,eSNP,Gpval from Geg, GSNP_eSNP_Gpval'
                        ' where Geg.GWAS = "' + GWAS + '"'
                        ' and Geg.eQTL = "' + eQTL + '"'
                        ' and Geg.gene = "' + gene + '"'
                        ' and GSNP_eSNP_Gpval.Geg_id = Geg.id;'  
                        )
        list_detail = self.exec_fetch_SQL(sql_template)
        return list_detail 
