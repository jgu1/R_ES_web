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
                
class Cluster(object):
    row_comb = None
    cols = None
    excluded = False 
    num_seeds = 1
    def __init__(self,row_comb,cols):
        self.row_comb = list(row_comb)
        self.cols = set(cols)
        self.excluded = False
        self.num_seeds = 1

class pair_manhattan(object):
    SNP_list = None
    def __init__(self,SNP_list):
        self.SNP_list = SNP_list

class gene_manhattan(object):
    SNP_list = None
    def __init__(self, SNP_list):
        self.SNP_list = SNP_list 


class SNP_manhattan(object):
    chrom = None
    pos   = None
    Gpval = None
    epval = None
    def __init__(self,chrom, pos, Gpval, epval):
        self.chrom  = chrom
        self.pos    = pos
        self.Gpval  = Gpval
        self.epval  = epval 

class Chrom_fields(object):
    Chrom_name_list = ['chr1' ,'chr2' ,'chr3' ,'chr4' ,'chr5' ,'chr6' ,'chr7' ,'chr8' ,'chr9' ,'chr10',
                       'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20',
                       'chr20','chr21','chr22','chrX' ,'chrY' ]
    Chrom_len_dict = {}
    Chrom_len_dict['chr1']  = 249250621 
    Chrom_len_dict['chr2']  = 243199373
    Chrom_len_dict['chr3']  = 198022430
    Chrom_len_dict['chr4']  = 191154276
    Chrom_len_dict['chr5']  = 180915260
    Chrom_len_dict['chr6']  = 171115067
    Chrom_len_dict['chr7']  = 159138663
    Chrom_len_dict['chr8']  = 146364022
    Chrom_len_dict['chr9']  = 141213431
    Chrom_len_dict['chr10'] = 135534747
    Chrom_len_dict['chr11'] = 135006516
    Chrom_len_dict['chr12'] = 133851895
    Chrom_len_dict['chr13'] = 115169878
    Chrom_len_dict['chr14'] = 107349540
    Chrom_len_dict['chr15'] = 102531392
    Chrom_len_dict['chr16'] = 90354753
    Chrom_len_dict['chr17'] = 81195210
    Chrom_len_dict['chr18'] = 78077248
    Chrom_len_dict['chr19'] = 59128983
    Chrom_len_dict['chr20'] = 63025520
    Chrom_len_dict['chr21'] = 48129895
    Chrom_len_dict['chr22'] = 51304566
    Chrom_len_dict['chrX']  = 155270560
    Chrom_len_dict['chrY']  = 59373566
    
    @staticmethod
    def calc_chr_start(chr_name):
        chr_start = 0
        chr_idx = Chrom_fields.Chrom_name_list.index(chr_name) 
        for idx in range(chr_idx):
            preceeding_chr_name = Chrom_fields.Chrom_name_list[idx]  # add the length of all preceeding chr 
            chr_start = chr_start + Chrom_fields.Chrom_len_dict[preceeding_chr_name]        
        return chr_start

class SNP_fields_raw_row_obj(object):
    _id             = None
    _gene_fields_id = None
    _block          = None
    _tag            = None
    _in_DB          = None 
    _proximity      = None
    _chrom          = None
    _eQTL_SNP       = None
    _PGC            = None
    _eQTL_location  = None
    _eQTL_MAF       = None
    _GWAS_SNP       = None
    _GWAS_location  = None
    _GWAS_MAF       = None
    _delta          = None
    _LD             = None
    _eQTL_pval      = None
    _GWAS_pval      = None

    def __init__(self,DB_row):
        self._id             = DB_row[0]
        self._gene_fields_id = DB_row[1]
        self._block          = DB_row[2]
        self._tag            = DB_row[3]
        self._in_DB          = DB_row[4] 
        self._proximity      = DB_row[5]
        self._chrom          = DB_row[6]
        self._eQTL_SNP       = DB_row[7]
        self._PGC            = DB_row[8]
        self._eQTL_location  = DB_row[9]
        self._eQTL_MAF       = DB_row[10]
        self._GWAS_SNP       = DB_row[11]
        self._GWAS_location  = DB_row[12]
        self._GWAS_MAF       = DB_row[13]
        self._delta          = DB_row[14]
        self._LD             = DB_row[15]
        self._eQTL_pval      = DB_row[16]
        self._GWAS_pval      = DB_row[17]

    def Err_msg(self,field_name):
        print 'ERROR converting "'+ field_name +'" for ' + str(self._id)

    def get_id(self):
        return self._id

    def get_gene_fields_id(self):
        return self._gene_fields_id

    def get_block(self):
        block = 0
        try:
            block = int(self._block)
        except ValueError:
            self.Err_msg('block')
        return block

    def get_tag(self):
        tag = (self._tag == 'Yes')
        return tag

    def get_in_DB(self):
        in_DB = (self._in_DB == 'Yes')
        return in_DB
    
    def get_proximity(self):
        proximity = self._proximity
        return proximity
      
    def get_chrom(self):
        chrom = self._chrom
        return chrom
        
    def get_eQTL_SNP(self):
        eQTL_SNP = self._eQTL_SNP
        return eQTL_SNP

    def get_PGC(self):
        PGC = None
        try:
            PGC = int(self._PGC)
        except ValueError:
            self.Err_msg('PGC')
        return PGC

    def get_eQTL_location(self):
        eQTL_location = None
        try:
            eQTL_location = int(self._eQTL_location) #sys.maxint = 9223372036854775807
        except ValueError:
            self.Err_msg('eQTL_location')       
        return eQTL_location
    
    def get_eQTL_MAF(self):
        eQTL_MAF = 0
        try:
            eQTL_MAF = float(self._eQTL_MAF)
        except ValueError:
            self.Err_msg('eQTL_MAF')
        return eQTL_MAF

    def get_GWAS_SNP(self):
        GWAS_SNP = self._GWAS_SNP
        if GWAS_SNP == 'NF':
            GWAS_SNP = None
        return GWAS_SNP
    
    def get_GWAS_location(self):
        GWAS_location = None
        try:
            GWAS_location = int(self._GWAS_location)
        except ValueError:
            self.Err_msg('GWAS_location') 
        return GWAS_location

    def get_GWAS_MAF(self):
        GWAS_MAF = 0
        try:
            GWAS_MAF = float(self._GWAS_MAF)
        except ValueError:
            self.Err_msg('GWAS_MAF')
        return GWAS_MAF

    def get_delta(self):
        delta = 0
        try:
            delta = int(self._delta)
        except ValueError:
            self.Err_msg('delta')
        return delta

    def get_LD(self):
        LD = 0
        try:
            LD = float(self._LD)
        except ValueError:
            self.Err_msg('LD')
        return LD

    def get_eQTL_pval(self):
        eQTL_pval = 1
        try:
            eQTL_pval = float(self._eQTL_pval)  # p = float('1.85e-300') works
        except ValueError:
            self.Err_msg('eQTL_pval')
        return eQTL_pval
   
    def get_GWAS_pval(self):
        GWAS_pval = 1
        try:
            GWAS_pval = float(self._GWAS_pval)
        except ValueError:
            self.Err_msg('GWAS_pval')
        return GWAS_pval
    
    def is_aligned(self):
        aligned = False
        if self.get_GWAS_SNP() is not None:
            aligned = True
        return aligned

    def is_tagged(self):
        return self.get_tag() 


    def gen_GWAS_tuple(self):
        GWAS_tuple = None
        
        GSNP_name       = self.get_GWAS_SNP()
        if GSNP_name is not None:   # return None if no aligned GWAS_SNP is found
            GSNP_chr        = self.get_chrom()
            GWAS_location   = self.get_GWAS_location()
            GSNP_pval       = self.get_GWAS_pval()
            aligned         = self.is_aligned()
            tagged          = self.is_tagged()

            GSNP_abs        = Chrom_fields.calc_chr_start(GSNP_chr) + GWAS_location
           
            GWAS_tuple = (GSNP_name,GSNP_chr,GSNP_abs,GSNP_pval,aligned,tagged)
        return GWAS_tuple

    def gen_eQTL_tuple(self,gene_name):
        eSNP_name           = self.get_eQTL_SNP()
        eSNP_chr            = self.get_chrom()
        eSNP_location       = self.get_eQTL_location()
        eSNP_pval           = self.get_eQTL_pval()
        aligned             = self.is_aligned()
        tagged              = self.is_tagged()
        GSNP_name           = self.get_GWAS_SNP()

        eSNP_abs            = Chrom_fields.calc_chr_start(eSNP_chr) + eSNP_location

        eQTL_tuple = (eSNP_name,eSNP_chr,eSNP_abs,eSNP_pval,aligned,tagged,gene_name,GSNP_name)
        return eQTL_tuple

