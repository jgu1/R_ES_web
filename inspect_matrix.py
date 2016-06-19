import pickle
import os,re,time
import numpy
import pdb
import itertools
from beans import Cluster
import pyRserve,math
import multiprocessing 
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
            if len(curr_potential_child_row_comb) == 1 or len(curr_potential_child_cols) == 1:
                curr_potential_child.excluded = True # if child has only one row or one col, exclude it
                continue

            if set(curr_potential_child_row_comb).issubset(set(curr_parent_row_comb)) and curr_potential_child_cols.issubset(curr_parent_cols):
                curr_potential_child.excluded = True # if child can be 'contained' in parent, exclude it
                curr_parent.num_seeds = curr_parent.num_seeds + 1
            # see how many seeds result in the same sub_cluster 
            #if set(curr_potential_child_row_comb) == set(curr_parent_row_comb) and curr_potential_child_cols== curr_parent_cols:
            #    curr_parent.num_seeds = curr_parent.num_seeds + 1
 
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
    for i_cluster in range(num_cluster):
        #pdb.set_trace()
        curr_cluster_row = rowCol['Row'][:,i_cluster]
        curr_cluster_col = rowCol['Col'][i_cluster]
        curr_cluster_row_index = [i_row for i_row ,x in enumerate(curr_cluster_row) if x]  # extract the index of true element
        #curr_cluster_row_name  = [disease_names[i_row] for i_row in curr_cluster_row_index]
        curr_cluster_col = [i_col for i_col,x in enumerate(curr_cluster_col) if x]
        
        curr_cluster_disease_names = [disease_names[i_col] for i_col in curr_cluster_col]
        curr_cluster = Cluster(curr_cluster_disease_names, curr_cluster_col)
        clusters.append(curr_cluster) 
        
    return clusters

def build_p_m_from_File():
    #IF_name = os.getcwd() + '/' + 'cluster_input_q0.5_log10p.txt' 
    IF_name = os.getcwd() + '/' + 'search_matrix.txt' 
    IF = open(IF_name,'r')
    diseases = IF.readline()
    diseases = re.split(r'\t',diseases)
    diseases = diseases[1:]
    num_disease = len(diseases)
    gene_pvals_d = dict()
    for line in IF:
        fields = re.split(r'\t',line)
        gene = fields[0]
        gene_pvals_d[gene] = fields[1:]
    num_gene = len(gene_pvals_d)
    matrix_dimension = (num_gene,num_disease)
    p_m = numpy.zeros(matrix_dimension)
    i_row = 0

    for gene, pvals_str in gene_pvals_d.iteritems():
        pvals_float = []
        for pval in pvals_str:
            curr_pval_float = 1e-8
            try:
                curr_pval_float = float(pval)
            except ValueError:
                a = 1
                #pdb.set_trace()    
            pvals_float.append(curr_pval_float)
        p_m[i_row,:] = pvals_float
        i_row = i_row + 1

    return p_m,diseases

def R_build_numpy_matrix_from_gene_p_qs(gene_p_qs,cutoff):
    #pdb.set_trace()
    num_row = len(gene_p_qs)
    keys = gene_p_qs.keys()
    num_col = len(gene_p_qs[keys[0]])
    print 'num_row = {} num_col = {}'.format(num_row,num_col)

    ndarr = numpy.zeros((num_row,num_col))

    #convert original matrix into a binary matrix, setting all elements passing cutoff to 1, other to 0
    for i_row in range(num_row):            
        curr_row = gene_p_qs[keys[i_row]]
        for i_col in range(num_col):
            curr_ele  = curr_row[i_col]
            curr_pval = float(curr_ele[3]) 
            if curr_pval>0 and curr_pval < cutoff:
                ndarr[i_row,i_col] = 1 
    
    num_pass = len(numpy.nonzero(ndarr)[0])    
    pass_rate = float(num_pass) /(num_row * num_col)
    return ndarr

#keep a row if it exceed per_cutoff std from mean or it contain more 1 than abs_cutoff
#converge if diff between prev_cols and curr_cols smaller than converge_epsilon or have iterated more than converge_depth times
#def manual_ISA(binary_mat,abs_cutoff, per_cutoff,converge_epsilon,converge_depth,seed0):
def manual_ISA(args):
    start_time = time.time()
    #pdb.set_trace()
    binary_mat       = args[0]
    abs_cutoff       = args[1]
    per_cutoff       = args[2]
    converge_epsilon = args[3]
    converge_depth   = args[4]
    seed0            = args[5] 
    #print '\nthis thread seed0 = '+ str(numpy.nonzero(seed0))

    num_row = binary_mat.shape[0]
    num_col = binary_mat.shape[1]
    if False:   #using threadPool, seed0 will be passed in
        conn = pyRserve.connect()
        conn.r('require("isa2")')
        seeds = conn.r('generate.seeds('+ str(num_col)+',count = 1)')
        num_seeds = seeds.shape[1]
        len_each_seed = seeds.shape[0]
        seed0 = seeds[:,0]

    prev_cols = seed0
    curr_depth = 0
    #print 'seed0 = ' + str(numpy.nonzero(seed0))
    while True:
        curr_rows = manual_ISA_filter_row(binary_mat,prev_cols,abs_cutoff,per_cutoff)
        curr_cols = manual_ISA_filter_col(binary_mat,curr_rows,abs_cutoff,per_cutoff)
        if not numpy.any(curr_cols):  #if converge to empty hole, terminate early
            print 'xxx ABORTION, ALL ZERO'
            break
        if converge(curr_cols,prev_cols,converge_epsilon):
            print '$$$ REAL CONVERGE'
            #pdb.set_trace()
            a = 1
            break
        elif curr_depth > converge_depth:
            print 'xxx TIME CONVERGE'
            #pdb.set_trace()
            a = 1
            break
        else:
            prev_cols = curr_cols   #iterate        
            curr_depth = curr_depth + 1
            
    #pdb.set_trace()
    #print 'this thread takes {} seconds'.format(time.time() - start_time)
    return curr_rows,curr_cols
            
     
#converge if      
def converge(curr_cols,prev_cols,converge_epsilon):
    diff_length = numpy.linalg.norm(numpy.subtract(curr_cols,prev_cols))
    sum_length  = numpy.linalg.norm(numpy.add(curr_cols,prev_cols))
    curr_epsilon = float(diff_length) / sum_length
    #print '###'
    #print 'curr_cols = ' + str(numpy.nonzero(curr_cols))
    #print 'prev_cols = ' + str(numpy.nonzero(prev_cols))
    #print 'curr_epsilon = ' + str(curr_epsilon)
    if curr_epsilon < converge_epsilon:
        return True
    else:
        return False
     
def manual_ISA_filter_row(binary_mat,col_rake,abs_cutoff,per_cutoff):
    num_row = binary_mat.shape[0]
    row_proj_result = numpy.zeros(num_row)
    for i_row in range(num_row):
        curr_row  = binary_mat[i_row,:]
        curr_proj = numpy.dot(curr_row,col_rake)
        row_proj_result[i_row] = curr_proj
    mean = numpy.mean(row_proj_result)
    std = numpy.std(row_proj_result)

    #keep a row if it exceed per_cutoff std from mean or it contain more 1 than abs_cutoff
    debug_count = 0
    filter_result = numpy.zeros(num_row)
    for i_row in range(num_row):
        curr_proj = row_proj_result[i_row]
        if curr_proj > mean + per_cutoff * std or curr_proj > abs_cutoff:
            filter_result[i_row] = 1
            debug_count = debug_count + 1
    return filter_result
         
def manual_ISA_filter_col(binary_mat,row_rake,abs_cutoff,per_cutoff):
    num_col = binary_mat.shape[1]
    col_proj_result = numpy.zeros(num_col)
    for i_col in range(num_col):
        curr_col  = binary_mat[:,i_col]
        curr_proj = numpy.dot(curr_col,row_rake)
        col_proj_result[i_col] = curr_proj
    mean = numpy.mean(col_proj_result)
    std = numpy.std(col_proj_result)

    #keep a row if it exceed per_cutoff std from mean or it contain more 1 than abs_cutoff
    debug_count = 0
    filter_result = numpy.zeros(num_col)
    for i_col in range(num_col):
        curr_proj = col_proj_result[i_col]
        if curr_proj > mean + per_cutoff * std or curr_proj > abs_cutoff:
            filter_result[i_col] = 1
            debug_count = debug_count + 1
    return filter_result

def manual_ISA_build_Cluster_objects(gene_p_qs,rows,cols):
    #pdb.set_trace()
    Cluster_row_comb = []
    disease_names = gene_p_qs.keys()
    for i_row,bool_i_row in enumerate(rows):
        if bool_i_row > 0:
            Cluster_row_comb.append(disease_names[i_row])
   
    idx_ndarr = numpy.nonzero(cols)
    Cluster_cols = idx_ndarr[0].tolist()
    #pdb.set_trace()
    return Cluster(Cluster_row_comb,Cluster_cols) 
        

def manual_ISA_gen_seeds(binary_mat,est_col_width,pre_exclude_gene_indices):
    num_row = binary_mat.shape[0]
    num_col = binary_mat.shape[1]
    num_seeds = int(num_col/est_col_width)   
    print 'in manual_ISA num_col = ' + str(num_col)
    print 'in manual_ISA num_seeds = ' + str(num_seeds)
    conn = pyRserve.connect()
    conn.r('require("isa2")')
    seeds_mat = conn.r('generate.seeds('+ str(num_col)+',count = '+str(num_seeds)+',sparsity='+str(est_col_width)+')')
    seeds_list = []
    num_seeds = seeds_mat.shape[1]

    #make a ndarray that having 1 at pre_excluded_indices
    pre_exclude_mask = numpy.zeros(num_col)
    for pre_exclude_gene_index in pre_exclude_gene_indices:
        pre_exclude_mask[pre_exclude_gene_index] = 1
 
    for i in range(num_seeds):
        curr_seed_vec = seeds_mat[:,i]
        #if curr_seed_vec contains pre_excluded gene_indices, remove this seed
        #the product is greater than 0 only when index matches
        if numpy.dot(curr_seed_vec,pre_exclude_mask) > 0: 
            continue
        seeds_list.append(curr_seed_vec)

    


    return seeds_list
    #num_seeds = seeds.shape[1]
    #len_each_seed = seeds.shape[0]
    #seed0 = seeds[:,0]

def manual_ISA_filter_sub_cluster(binary_mat, rows, cols, row_cutoff,col_cutoff):
    num_row = numpy.count_nonzero(rows)
    num_col = numpy.count_nonzero(cols)


    nonzero_col_indices = numpy.nonzero(cols)[0].tolist() #numpy.nonzero(cols) return a 2-elements tuple
    for i_col in nonzero_col_indices:
        curr_col = binary_mat[:,i_col]
        curr_sum = numpy.sum(curr_col)
        if curr_sum < num_row * col_cutoff:
            cols[i_col] = 0

def R_discover_sub_clusters(gene_p_qs,abs_cutoff,per_cutoff,converge_epsilon,converge_depth,est_col_width,filter_ratio):
    binary_mat = R_build_numpy_matrix_from_gene_p_qs(gene_p_qs,1E-3)
    #abs_cutoff = 3
    #per_cutoff = 0.5
    #converge_epsilon = 0.1
    #converge_depth = 100


    #pre_exclude_gene_names = ['DND1','LRRC37A4','MAPK8IP1','MAPT','ZNF285','CLDN23']
    pre_exclude_gene_names = []
    pre_exclude_gene_indices = get_pre_exclude_gene_idx(gene_p_qs,pre_exclude_gene_names)

    seeds = manual_ISA_gen_seeds(binary_mat,est_col_width,pre_exclude_gene_indices)
    output_seed_to_txt(seeds,gene_p_qs) 
    manual_ISA_args = []
    for seed in seeds:
        curr_arg = (binary_mat,abs_cutoff,per_cutoff,converge_epsilon,converge_depth,seed)
        manual_ISA_args.append(curr_arg)


    threadPool = multiprocessing.Pool(5)
    start_time = time.time()
    sub_clusters_rows_cols = threadPool.map(manual_ISA,manual_ISA_args)        
    print 'all threads take {} seconds'.format(time.time() - start_time)

    sub_clusters = []
    for rows_cols in sub_clusters_rows_cols:
        rows = rows_cols[0]
        cols = rows_cols[1]
        manual_ISA_filter_sub_cluster(binary_mat, rows, cols, filter_ratio, filter_ratio)
        if not numpy.any(rows) or not numpy.any(cols):
            continue

        #print '\n'
        #print 'len(rows) = ' + str(numpy.count_nonzero(rows))
        #print 'len(cols) = ' + str(numpy.count_nonzero(cols))



        curr_sub_cluster = manual_ISA_build_Cluster_objects(gene_p_qs,rows,cols)
        sub_clusters.append(curr_sub_cluster)


    #rows,cols = manual_ISA(binary_mat,abs_cutoff, per_cutoff,converge_epsilon,converge_depth)
    #one_Cluster = manual_IRA_build_Cluster_objects(gene_p_qs,rows,cols) 
    #return [one_Cluster]

    filtered_sub_clusters = filter_out_child_sub_clusters(sub_clusters)  
    ##merged_sub_clusters = merge_sub_clusters_with_same_cols(filtered_sub_clusters)
 
    return filtered_sub_clusters

    if False:
        start_time = time.time() 
        p_m = R_build_matrix(gene_p_qs)
        conn = pyRserve.connect()
        conn.r('require("biclust")')
        R_args = {
            'x':p_m,
            'method':'BCPlaid',
            'cluster':'b',
            #'fit.model':'y~m+a+b',
            'background':False,
            'row.release':0.7,
            'col.release':0.7,
            #'shuffle':3,
            'shuffle':3,
            'back.fit':0,
            'max.layers':20,
            'iter.startup':5,
            'iter.layer':10,
            'verbose':True,
        }
        result = conn.r.biclust(**R_args)
        attr = result.lexeme.attr
        disease_names = gene_p_qs.keys() #FIXME temporarily comment off for testing p_m
        clusters = R_parse_cluster_result(attr,disease_names)    

        print 'found ' + str(len(clusters)) + ' clusters\n'    
        clusters = R_filter_clusters(clusters,gene_p_qs,row_percent,row_cutoff,col_percent,col_cutoff)
        print("sub_clustering took --- %s seconds ---" % (time.time() - start_time))
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
            p_val_str = gene_p_qs[disease][i_gene][3]
            p_val = -1
            try:
                pval = float(p_val_str)
                pval =  -numpy.log10(pval)
            except ValueError:
                pval = 1e-8
            OF.write(str(pval)) 

def get_pre_exclude_gene_idx(gene_p_qs,pre_exclude_gene_names):
    pre_exclude_gene_indices = set()
    first_gene_tuple_list = gene_p_qs.values()[0]
    for gene_idx,gene_tuple in enumerate(first_gene_tuple_list):
        curr_gene_name = gene_tuple[2]
        if curr_gene_name in pre_exclude_gene_names:
            pre_exclude_gene_indices.add(gene_idx)
    return pre_exclude_gene_indices
            
def output_seed_to_txt(seeds,gene_p_qs):   

    SEED_TXT = open('seed_genes.txt','w+')

    first_gene_tuple_list = gene_p_qs.values()[0]
    all_gene_names = [x[2] for x in first_gene_tuple_list]
    all_gene_names.sort()
    for seed in seeds:

        # print genes in current seed
        curr_seed_genes = []
        seed_length = seed.shape[0]
        for i in range(seed_length):
            if seed[i] > 0:
                curr_seed_genes.append(all_gene_names[i])
        SEED_TXT.write( '\ncurr_seed_genes = ' + str(curr_seed_genes) )
        # print genes in current seed
  

if __name__=='__main__':
    pickle_filename='ES_Sherlock_dump.pickle'
    Main(os.getcwd() + '/' + pickle_filename)
