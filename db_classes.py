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
