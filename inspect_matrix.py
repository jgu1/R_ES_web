import pickle
import os
import numpy
import pdb
import itertools
from beans import Cluster
def Main(pickle_filename):
    matrix = pickle.load(open(pickle_filename))
    pdb.set_trace()
    pass_dict = build_pass_dict(matrix)
    row_combs = build_row_combinations(matrix)
    sub_clusters = detect_sub_clusters(row_combs,pass_dict)
    pdb.set_trace()
  
    filtered_sub_clusters = filter_out_child_sub_clusters(sub_clusters)  
    pdb.set_trace()
    a = 1

def discover_sub_clusters(matrix):
    pass_dict = build_pass_dict(matrix)
    row_combs = build_row_combinations(matrix)
    sub_clusters = detect_sub_clusters(row_combs,pass_dict,matrix)
    filtered_sub_clusters = filter_out_child_sub_clusters(sub_clusters)  
    pdb.set_trace()
    merged_sub_clusters = merge_sub_clusters_with_same_cols(filtered_sub_clusters)
    return merged_sub_clusters

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
    sub_clusters = sorted(sub_clusters, key=lambda x: len(x.row_comb), reverse=True)
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

if __name__=='__main__':
    pickle_filename='ES_Sherlock_dump.pickle'
    Main(os.getcwd() + '/' + pickle_filename)
