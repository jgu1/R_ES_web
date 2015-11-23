import os
import pdb
class Gene_p_q(object):
    GWAS=None
    eQTL=''
    gene=None
    pval=''
    qval=''
    _allGWAS=[]
    _alleQTL=[]

    def read_all_lines(self,GWAS_or_eQTL):
        filename = None
        if GWAS_or_eQTL == 'GWAS':
            filename = os.getcwd() + '/allGWAS'
        elif GWAS_or_eQTL == 'eQTL':
            filename = os.getcwd() + '/alleQTL'
               
        with open(filename,'r') as allGWAS_file:
            GWAS_lines = allGWAS_file.readlines()
            for line in GWAS_lines:
                line = line.strip()
                if not not line:
                    if GWAS_or_eQTL == 'GWAS':
                        self._allGWAS.append(line)
                    elif GWAS_or_eQTL == 'eQTL':
                        self._alleQTL.append(line)             

    def parse_sherlock_output_pairname(self,sherlock_output_pairname):
        pairname = ''
        if sherlock_output_pairname.endswith('_gene.txt'):
            pairname = sherlock_output_pairname[:-9]
        for splitPos in range(1,len(pairname)):
            if pairname[:splitPos] in self._allGWAS and pairname [splitPos+1:] in self._alleQTL:
                break
        GWAS = pairname[:splitPos]
        eQTL = pairname[splitPos+1:]
        return GWAS,eQTL        


    def __init__(self,sherlock_output_pairname,gene,pval,qval):
        self._allGWAS = []
        self._alleQTL = []    
        self.read_all_lines('GWAS')
        self.read_all_lines('eQTL')

        GWAS,eQTL = self.parse_sherlock_output_pairname(sherlock_output_pairname)
        self.GWAS = GWAS
        self.eQTL = eQTL

        self.gene = gene
        self.pval = pval
        self.qval = qval
                

class Gene_SNPs(object):
    GWAS=None
    eQTL=None
    gene=None
    GSNP_eSNP_Gpvals = []

    def __init__(self,GWAS,eQTL,gene,GSNP_eSNP_Gpvals):
        self.GWAS = GWAS
        self.eQTL = eQTL
        self.gene = gene
        self.GSNP_eSNP_Gpvals = GSNP_eSNP_Gpvals
                














