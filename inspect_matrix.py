import pickle
import os
import numpy
import pdb
import itertools
from beans import Cluster
import pyRserve,math
OF_TXT_NAME = os.getcwd() + '/search_matrix.txt'

def Main(pickle_filename):
#    matrix = pickle.load(open(pickle_filename))
#    pdb.set_trace()
#    pass_dict = build_pass_dict(matrix)
#    row_combs = build_row_combinations(matrix)
#    sub_clusters = detect_sub_clusters(row_combs,pass_dict)
  
#    filtered_sub_clusters = filter_out_child_sub_clusters(sub_clusters)  
#    merged_sub_clusters = merge_sub_clusters_with_same_cols(filtered_sub_clusters)

    gene_p_qs = pickle.load(open(pickle_filename))
    pdb.set_trace()
    clusters = R_discover_sub_clusters(gene_p_qs)
    pdb.set_trace()
    a = 1


def discover_sub_clusters(matrix):
    pass_dict = build_pass_dict(matrix)
    row_combs = build_row_combinations(matrix)
    sub_clusters = detect_sub_clusters(row_combs,pass_dict,matrix)
    #return sub_clusters
    filtered_sub_clusters = filter_out_child_sub_clusters(sub_clusters)  
    #return filtered_sub_clusters
    merged_sub_clusters = merge_sub_clusters_with_same_cols(filtered_sub_clusters)
    #return merged_sub_clusters
    sub_clusters = filter_out_child_sub_clusters(merged_sub_clusters)  
    return sub_clusters

def merge_sub_clusters_with_same_cols(sub_clusters):
    dict_cols_2_row_combs={}
    for sub_cluster in sub_clusters:
        cols_tuple = tuple(sub_cluster.cols)
        curr_row_comb_set = set(sub_cluster.row_comb)
        if cols_tuple in dict_cols_2_row_combs:
            dict_cols_2_row_combs[cols_tuple] = dict_cols_2_row_combs[cols_tuple].union(curr_row_comb_set)
        else:
            dict_cols_2_row_combs[cols_tuple] = curr_row_comb_set
   
    merged_sub_clusters = []

    for key, value in dict_cols_2_row_combs.iteritems():
        cols = set(key)
        row_comb = list(value)
        row_comb.sort()
        cluster = Cluster(row_comb, cols)
        merged_sub_clusters.append(cluster)
    
    return merged_sub_clusters

def filter_out_child_sub_clusters(sub_clusters):
    # sort by number of rows, if same number of rows, sort by number of cols
    sub_clusters = sorted(sub_clusters, key=lambda x: (len(x.row_comb),len(x.cols)), reverse=True)
    for i in xrange(len(sub_clusters)-1):   # doesn't need to consider the last one as parent
        curr_parent = sub_clusters[i]
        if curr_parent.excluded == True:
            continue
        curr_parent_row_comb = curr_parent.row_comb
        curr_parent_cols = curr_parent.cols
        for j in xrange(i + 1, len(sub_clusters)):  #looking for child from the next cluster
            curr_potential_child = sub_clusters[j]
            curr_potential_child_row_comb =  curr_potential_child.row_comb
            curr_potential_child_cols = curr_potential_child.cols
            if set(curr_potential_child_row_comb).issubset(set(curr_parent_row_comb)) and curr_potential_child_cols.issubset(curr_parent_cols):
                curr_potential_child.excluded = True # if child can be 'contained' in parent, exclude it
    
    filtered_sub_clusters = []
    for cluster in sub_clusters:
        if cluster.excluded == False:
            filtered_sub_clusters.append(cluster)
    return filtered_sub_clusters

def detect_sub_clusters(row_combs,pass_dict,matrix):
    sub_clusters = []
    # scan through all row combinations
    for row_comb in row_combs:
        # for a column to be included, it must have over 60% percent of rows pass the threshold of pval 0.001
        cols = find_majority_cols(row_comb,pass_dict)
        if not not cols :
            
            row_comb_containing_at_least_one_pass_col = screen_rows_given_cols(row_comb, cols, pass_dict)
            cluster = Cluster(row_comb_containing_at_least_one_pass_col,cols) # there can be duplicating clusters, but that's fine because duplication will be removed by filter_out_child_sub_clusters
            '''
            cluster = Cluster(row_comb, cols) # if not screening rows, many similar cluster will be spawn because every row can become 'drive-by' row along with 'all-red' rows
            '''
            sub_clusters.append(cluster)
    return sub_clusters
   
def screen_rows_given_cols(row_comb, cols, pass_dict):
    row_comb_containing_at_least_one_pass_col = []
    for row in row_comb:
        pass_cols_for_curr_row = pass_dict[row]
        if len(pass_cols_for_curr_row.intersection(cols)) > 0:
            row_comb_containing_at_least_one_pass_col.append(row)
        '''
        else:
            pdb.set_trace()
            a = 1
        '''
    return row_comb_containing_at_least_one_pass_col
 
def find_comm_cols(row_comb,pass_dict):
    all_sets = pass_dict.values()
    comm_cols = set().union(*all_sets)
    for row in row_comb:
        curr_row_passes = pass_dict[row]
        comm_cols = comm_cols.intersection(curr_row_passes)
        if len(comm_cols) == 0:
            break
    if len(comm_cols) > 0:
       return comm_cols
    else:
        return None

def find_majority_cols(row_comb,pass_dict):
    cutoff = 0.6 # a col needs to have pval greater than 3 in 60% rows to be included
    all_sets = pass_dict.values()
    all_cols = set().union(*all_sets)
    num_rows = len(row_comb)
    majority_cols = set()
    for col in all_cols:
        pass_count = 0
        for row in row_comb:
            if col in pass_dict[row]:
               pass_count = pass_count + 1 
        if float(pass_count)/num_rows > cutoff:
            majority_cols.add(col)
    return majority_cols 

# return all row combinations, combination length starts from 2, to the whole length 
def build_row_combinations(matrix):
    keys = matrix.keys()
    combs = []
    for comb_length in range(2,len(keys) + 1):
        els = [list(x) for x in itertools.combinations(keys, comb_length)]
        combs = combs + els
    return combs
 
# screen each row, return idx 
def build_pass_dict(matrix):
    pass_dict={}
    for pair in matrix:
        pair_list = matrix[pair]
        pair_idx_list=[]
        for i in range(len(pair_list)):
            curr_tuple = pair_list[i]
            if pass_cutoff(curr_tuple):
                pair_idx_list.append(i)
        #pair_idx_tuple = tuple(pair_idx_list)
        pair_idx_set = set(pair_idx_list)
        pass_dict[pair] = pair_idx_set
    return pass_dict 

def pass_cutoff(tup):
    cutoff = 3
    pval = float(tup[3])
    if pval < 0:
        return False
    else:
        return -numpy.log10(pval) > cutoff

def pass_sparse_check(row_comb,gene_idx_set,matrix):
    cutoff = 0.75
    num_row = len(row_comb)
    num_gene = len(gene_idx_set)
    num_tot = num_row * num_gene
    num_pass = 0
    for row in row_comb:
        curr_row = matrix[row]
        for gene_idx in gene_idx_set:
            curr_tuple = curr_row[gene_idx]
            if pass_cutoff(curr_tuple):
                num_pass = num_pass + 1
    if float(num_pass)/num_tot >= cutoff:
        return True
    else:
        return False 

def R_build_matrix(gene_p_qs):
    disease_names = gene_p_qs.keys()
    disease_names.sort() 
    num_disease = len(gene_p_qs)
    num_gene = -1
    #loop one iteration to get number of genes
    for key,value in gene_p_qs.iteritems():
        num_gene = len(value)
        break
    matrix_dimension = (num_disease,num_gene)
    p_m = numpy.zeros(matrix_dimension)
    for i_disease, disease_name in enumerate(disease_names):
        curr_genes = gene_p_qs[disease_name]
        for i_gene,gene in enumerate(curr_genes):
            curr_pval_str = gene[3]
            curr_pval_float = 1e-8
            try:
                curr_pval_float = float(curr_pval_str)
            except ValueError:
                a = 1
            p_m[i_disease,i_gene] = curr_pval_float
    return p_m

    
def R_filter_clusters(clusters,gene_p_qs,row_percent,row_cutoff,col_percent,col_cutoff):
    new_clusters = []
    for curr_cluster in clusters:
        new_row_combs = []
        for pair_name in curr_cluster.row_comb:
            if R_filter_cluster_one_row(gene_p_qs,pair_name,curr_cluster.cols,row_percent,row_cutoff):
                new_row_combs.append(pair_name)
        
        new_cols = R_filter_cluster_cols(gene_p_qs,new_row_combs, curr_cluster.cols, col_percent,col_cutoff)

        curr_new_cluster = Cluster(new_row_combs,new_cols)
        new_clusters.append(curr_new_cluster)
    return new_clusters

#for the given row, if row_cutoff % column pass threshold, keep this row
def R_filter_cluster_one_row(gene_p_qs,pair_name, cols, row_percent, row_cutoff):
    all_genes = gene_p_qs[pair_name]
    num_cols = len(cols)
    num_pass = 0
    for i_gene in cols:
        curr_gene = all_genes[i_gene]
        curr_pval = float(curr_gene[3])
        if curr_pval < 0 or curr_pval > row_cutoff:
            continue
        else:
            num_pass = num_pass + 1
    if num_pass >= num_cols * row_percent:
        return True
    else:
        return False

def R_filter_cluster_cols(gene_p_qs,new_row_combs, cols, col_percent,col_cutoff):
    new_cols = []
    num_rows = len(new_row_combs)
    for i_gene in cols:
        num_pass = 0
        for pair_name in new_row_combs: 
            curr_row = gene_p_qs[pair_name]
            curr_gene = curr_row[i_gene]
            #if (curr_gene[2] == 'ABCA8'):
            #    pdb.set_trace()
            curr_pval = float(curr_gene[3])
            if curr_pval < 0 or curr_pval > col_cutoff:
                continue
            else:
                num_pass = num_pass + 1
        if num_pass >= num_rows * col_percent:            
            new_cols.append(i_gene)
    return new_cols

def R_parse_cluster_result(attr,disease_names):
    rowCol = dict()
    num_cluster = -1
    clusters = []
    for elem in attr:
        if elem[0] == 'RowxNumber':
            rowIndex = elem[1]
            rowCol['Row'] = rowIndex
        if elem[0] == 'NumberxCol':
            colIndex = elem[1]
            rowCol['Col'] = colIndex
            num_cluster = len(colIndex)
    for i in range(num_cluster):
        curr_cluster_row = rowCol['Row'][:,i]
        curr_cluster_col = rowCol['Col'][i]
        curr_cluster_row_index = [i for i,x in enumerate(curr_cluster_row) if x]  # extract the index of true element
        curr_cluster_row_name  = [disease_names[i] for i in curr_cluster_row_index]
        curr_cluster_col = [i for i,x in enumerate(curr_cluster_col) if x]
        
        curr_cluster = Cluster(curr_cluster_row_name, curr_cluster_col)
        clusters.append(curr_cluster) 
        
    return clusters

def R_discover_sub_clusters(gene_p_qs,row_percent,row_cutoff,col_percent,col_cutoff):
    p_m = R_build_matrix(gene_p_qs)
    conn = pyRserve.connect()
    conn.r('require("biclust")')
    if False:
        R_args = {
            'x':p_m,
            'method':'BCPlaid',
            'cluster':'b',
            'background':False,
            'row.release':0.7,
            'col.release':0.7,
            'shuffle':19,
            'back.fit':0,
            'max.layers':20,
            'iter.startup':5,
            'iter.layer':10,
            'verbose':True,
        }
        result = conn.r.biclust(**R_args)
    result = conn.r.biclust(p_m, method = "BCPlaid",cluster = 'b',background = False, shuffle = 19, verbose = True)
    attr = result.lexeme.attr
    disease_names = gene_p_qs.keys()
    clusters = R_parse_cluster_result(attr,disease_names)    
    
    #row_percent = 0.1
    #row_cutoff = 1E-1
    #col_percent = 0.3
    #col_cutoff = 1E-2 
    clusters = R_filter_clusters(clusters,gene_p_qs,row_percent,row_cutoff,col_percent,col_cutoff)
    return clusters

def output_matrix_to_txt(ret):
    genes = ret['filtered_gene_names']
    diseases = ret['sorted_pair_names']
    gene_p_qs = ret['gene_p_qs']
    OF = open(OF_TXT_NAME,'w+')
    #write the header
    OF.write('gwas')
    for d in diseases:
        OF.write('\t')
        OF.write(d)
    #write the header
   
    
    for i_gene,gene in enumerate(genes):    #write each row
        OF.write('\n')
        OF.write(gene)
        for i_d,disease in enumerate(diseases):  #write each column
            OF.write('\t')
            OF.write(gene_p_qs[disease][i_gene][3]) 

if __name__=='__main__':
    pickle_filename='ES_Sherlock_dump.pickle'
    Main(os.getcwd() + '/' + pickle_filename)
