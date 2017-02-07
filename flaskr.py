#all the imports

import pdb
import MySQLdb
import os
from flask import Flask, request, session, g, redirect, url_for, \
     abort, render_template, flash, jsonify, Response
from contextlib import closing
from time import sleep
from flask.ext.paginate import Pagination
import re
import socket
from datetime import datetime,timedelta
from db_classes import DAO
import json
import pickle
from inspect_matrix import R_discover_sub_clusters_ISA, R_discover_sub_clusters_PLAID, discover_sub_clusters,output_matrix_to_txt
import time
# configuration
DEBUG = True
SECRET_KEY = 'development key'
eQTL_name_all = 'all'
eQTL_names = ['Dixon_07','Duan_08','Liang_2012','Muther_12','Myers_07','Schadt_08','Wright_14_pruned_e6_2','Zeller_10','merged_pickle','GTExV6_LD85']
#disease_GWAS_dict = pickle.load(open("disease_GWAS_dict.pickle","rb"))
disease_GWAS_dict = DAO(None,None).build_disease_GWAS_dict()
GENE_P_Q_PER_PAGE=30
page=1

app = Flask(__name__)
app.config.from_object(__name__)
def connect_db():

    if 'web_disease_list' not in session or 'web_eQTL_list' not in session:
        return DAO(None,None)
    web_disease_list = session['web_disease_list']
    web_eQTL_list = session['web_eQTL_list']
    return DAO(web_disease_list, web_eQTL_list)

@app.before_request
def before_request():
    g.dao = connect_db()

@app.teardown_request
def teardown_request(exception):
    dao = getattr(g, 'dao', None)
    if dao is not None:
        dao.db.close()

@app.route('/sub_clusters')
@app.route("/sub_clusters/<string:alg><string:filter_ratio><string:binarize_cutoff><string:consider_all_genes_in_database>")
@app.route("/sub_clusters/<string:alg><string:binarize_cutoff><string:consider_all_genes_in_database>")
def sub_clusters():
    binarize_cutoff  = request.args.get('binarize_cutoff')
    consider_all_genes_in_database = request.args.get('consider_all_genes_in_database','')
    if consider_all_genes_in_database == 'true':
        consider_all_genes_in_database = True
    else:
        consider_all_genes_in_database = False
 
    alg = request.args.get('alg','')
    if alg == 'ISA':
        
        abs_cutoff       = request.args.get('abs_cutoff', '')
        per_cutoff       = request.args.get('per_cutoff','')
        converge_epsilon = request.args.get('converge_epsilon','')
        converge_depth   = request.args.get('converge_depth','')
        est_col_width    = request.args.get('est_col_width','')
        filter_ratio     = request.args.get('filter_ratio','')
        if abs_cutoff == '':
            abs_cutoff = 3
        if per_cutoff == '':
            per_cutoff  = 0.4
        if converge_epsilon == '':
            converge_epsilon = 0.1
        if converge_depth == '':
            converge_depth  = 100
        if est_col_width == '':
            est_col_width = 20
        if binarize_cutoff == '':
            binarize_cutoff = 0.003157#1E-2.5
        if filter_ratio == '':
            filter_ratio = 0.5
   
        gene_p_qs,filtered_gene_names,gene_descriptions = fetch_and_build_matrix(consider_all_genes_in_database)
        sub_clusters = R_discover_sub_clusters_ISA(gene_p_qs,float(abs_cutoff),float(per_cutoff),float(converge_epsilon),float(converge_depth),float(est_col_width),float(filter_ratio),float(binarize_cutoff))
   
    
    elif alg == 'PLAID':
        if binarize_cutoff == '':
            binarize_cutoff = 0.0001
        gene_p_qs,filtered_gene_names,gene_descriptions = fetch_and_build_matrix(consider_all_genes_in_database)
        sub_clusters = R_discover_sub_clusters_PLAID(gene_p_qs,float(binarize_cutoff))

    dao = getattr(g, 'dao', None)
    sub_clusters = dao.remove_sub_clusters_from_same_disease(sub_clusters)
 
    serisables = []
    for sub_cluster in sub_clusters:
        row_comb = sub_cluster.row_comb
        cols = list(sub_cluster.cols)
        cols.sort()
        num_seeds = sub_cluster.num_seeds
        curr_serisable = [row_comb,cols,num_seeds]
        serisables.append(curr_serisable)


    ret = {}
    ret['gene_p_qs'] = gene_p_qs
    ret['gene_descriptions'] = gene_descriptions
    ret['serisables'] = serisables

    resp = Response(json.dumps(ret), status=200, mimetype='application/json')
    return resp   
 

@app.route('/detail')
@app.route("/data/<string:gene>")
@app.route("/data/<string:gene><string:pairNames>")
def detail():
    gene = request.args.get('gene', 'empty')
    pairNames = request.args.get('pairNames','empty')
    web_disease_list = session['web_disease_list']
    web_eQTL_list = session['web_eQTL_list']
    dao = getattr(g, 'dao', None)
    pair_SNP_dict_all,all_SNPs_list = dao.fetch_pair_SNP(web_disease_list,web_eQTL_list,gene)
 
    pair_SNP_dict = {}
    if pairNames != 'empty':
        for pair in pairNames.split(','):
            if pair in pair_SNP_dict_all:
                pair_SNP_dict[pair] = pair_SNP_dict_all[pair]
            else:   # this gene may not be in the current pair. 
                continue
    else:
        pair_SNP_dict = pair_SNP_dict_all

    #pair_SNP_dict = pair_SNP_dict_all
    ret = {}
    ret['gene'] = gene
    ret['pair_SNP_dict'] = pair_SNP_dict
    ret['all_SNPs_list'] = all_SNPs_list
    return jsonify(ret)


def Manhattan_gen_pair_SNP_dict(web_disease_list,web_eQTL_list,pairNames,gene):
    dao = getattr(g, 'dao', None)
    pair_SNP_dict_all,all_SNPs_list = dao.fetch_pair_SNP(web_disease_list,web_eQTL_list,gene)
    pair_SNP_dict = {}
    if pairNames != 'empty':
        for pair in pairNames.split(','):
            if pair in pair_SNP_dict_all:
                pair_SNP_dict[pair] = pair_SNP_dict_all[pair]
            else:   #this gene may not appear in this pair, in this case no SNP can be ever found
                continue
    else:
        pair_SNP_dict = pair_SNP_dict_all
    return pair_SNP_dict

def Manhattan_add_pair_SNP_dict(pair_SNP_dict_with_location_all_genes,pair_SNP_dict_with_location_curr_gene):
    pair_names_curr_gene = pair_SNP_dict_with_location_curr_gene.keys()
    
    for pair_name in pair_names_curr_gene:
        
        if pair_name not in pair_SNP_dict_with_location_all_genes:
            pair_SNP_dict_with_location_all_genes[pair_name] = pair_SNP_dict_with_location_curr_gene[pair_name] 
        else:
            existing_SNP_dict = pair_SNP_dict_with_location_all_genes[pair_name]
            pair_SNP_dict_with_location_all_genes[pair_name] = existing_SNP_dict + pair_SNP_dict_with_location_curr_gene[pair_name]
            
    return pair_SNP_dict_with_location_all_genes

# the input 'GWAS_SNPlist_dict_for_curr_gene' is a dict, has key as 'GWAS_eQTL', value as a list of SNPs for current GWAS from curr_gene
# this function adds all GWAS SNPs from current gene into the SNPlist of certain pair from ALL genes 
def Manhattan_add_GWAS_SNPlist_dict(GWAS_SNPlist_dict,GWAS_SNPlist_dict_for_curr_gene):
    pair_names_curr_gene = GWAS_SNPlist_dict_for_curr_gene.keys()
    
    for pair_name in pair_names_curr_gene:
        # each pair has a list of SNPs
        if pair_name not in GWAS_SNPlist_dict:
            GWAS_SNPlist_dict[pair_name] = GWAS_SNPlist_dict_for_curr_gene[pair_name]
        else:
            existing_SNPlist = GWAS_SNPlist_dict[pair_name]
            GWAS_SNPlist_dict[pair_name] = existing_SNPlist + GWAS_SNPlist_dict_for_curr_gene[pair_name]

    return GWAS_SNPlist_dict

# the input 'eQTL_SNPlist_dict_for_curr_gene' is a dict, has key as 'GWAS_eQTL', value as a list of SNPs from current eQTL from curr_gene
# inside the 'container' data structure eQTL_gene_SNPlist_dict, there are 2 levels of dictionary_ing
# first level has key as 'GWAS_eQTL', value as a dict
    # the nested dict has key as 'gene_name', value as a list of SNPS belonging to specific gene and specific eQTL

# the input eQTL_SNPlist_dict_for_curr_gene has only 1 level of dictionary_ing
# key is 'GWAS_eQTL'
# value is a list of eQTL SNPs for this particular gene and particular eQTL

# this function first untangles the first level of dict, get a 'gene_SNPlist_dict_for_curr_pair'
    # this is a dict with key as 'gene_name',value as a list of SNPS belonging to specific gene and specific eQTL
    # if no genes have been added to the curr pair, a empty dict is created. This happens only once for each GWAS_eQTL pair
# then this function adds SNP_list for current gene and eQTL into 'gene_SNPlist_dict_for_curr_pair' 
def Manhattan_add_eQTL_gene_SNPlist_dict(eQTL_gene_SNPlist_dict,eQTL_SNPlist_dict_for_curr_gene,curr_gene_name):
    pair_names_curr_gene = eQTL_SNPlist_dict_for_curr_gene.keys()

    for pair_name in pair_names_curr_gene:
        #each pair has a dict{gene:SNPlist}
        gene_SNPlist_dict_for_curr_pair = {}
        try:
            gene_SNPlist_dict_for_curr_pair = eQTL_gene_SNPlist_dict[pair_name]

        except KeyError:# this happens only once for each GWAS_eQTL pair
            eQTL_gene_SNPlist_dict[pair_name] = gene_SNPlist_dict_for_curr_pair 
            # at this moment this is simply an empty dict             
            # use dict name instead of {} because we will add stuff into the empty dict in the immediate following code 
        eSNP_list_for_particular_gene_and_particular_eQTL = eQTL_SNPlist_dict_for_curr_gene[pair_name]
        gene_SNPlist_dict_for_curr_pair[curr_gene_name] = eSNP_list_for_particular_gene_and_particular_eQTL
           
    return eQTL_gene_SNPlist_dict

def Manhattan_debug(pair_SNP_dict_with_location_all_genes):
    location_pval_chrom_SNPlist_dict = {}
    eQTL_SNPlist_dict = {}
    for pair_name,SNP_tuple_list_all in pair_SNP_dict_with_location_all_genes.iteritems():
        curr_location_pval_chrom_SNPlist = []
        curr_eQTL_SNPlist = []
        for SNP_tuple in SNP_tuple_list_all:
            try:
                GSNP_name       = SNP_tuple[0]
                eSNP_name       = SNP_tuple[1]
                if GSNP_name == 'dummy' or eSNP_name == 'dummy':    # skip any dummy tuples
                    continue
                GSNP_pval       = SNP_tuple[2]
                eSNP_pval       = SNP_tuple[3]
                GSNP_chrom_abs  = SNP_tuple[4]
                GSNP_chrom      = GSNP_chrom_abs[0]
                GSNP_abs        = GSNP_chrom_abs[1]
                eSNP_chrom_abs  = SNP_tuple[5]
                eSNP_chrom      = eSNP_chrom_abs[0]
                eSNP_abs        = eSNP_chrom_abs[1]
                gene            = SNP_tuple[6]
            except TypeError:
                pdb.set_trace()
                a = 1
     
            curr_location_pval_chrom_SNPlist_tuple = (GSNP_name,GSNP_abs,GSNP_pval,gene,GSNP_chrom)
            curr_eQTL_SNPlist_tuple = (eSNP_name,eSNP_abs,eSNP_pval,gene,eSNP_chrom) 
             
            curr_location_pval_chrom_SNPlist.append(curr_location_pval_chrom_SNPlist_tuple)
            curr_eQTL_SNPlist.append(curr_eQTL_SNPlist_tuple)

        location_pval_chrom_SNPlist_dict[pair_name] = curr_location_pval_chrom_SNPlist
        eQTL_SNPlist_dict[pair_name] = curr_eQTL_SNPlist

    return location_pval_chrom_SNPlist_dict,eQTL_SNPlist_dict

@app.route('/Manhattan')
@app.route("/Manhattan/<string:geneNames><string:pairNames><string:GSNP_cutoff>")
def Manhattan():
    start_time_Manhattan = time.time()
    geneNames = request.args.get('geneNames', 'empty')
    pairNames = request.args.get('pairNames','empty')
    GSNP_cutoff_str = request.args.get('GSNP_cutoff','empty')

    GSNP_cutoff = None
    try:
        GSNP_cutoff = float(GSNP_cutoff_str)
    except ValueError:
        a = 1  # placeholder 


    web_disease_list = session['web_disease_list']
    web_eQTL_list = session['web_eQTL_list']

    genes = geneNames.split(',')
    genes = filter(None, genes) #remove empty gene, empty gene is unchecked in the web_interface   
    genes.sort()
 
    dao = getattr(g, 'dao', None)
    location_pval_chrom_SNPlist_dict = {}
    Manhattan_pairNames = set()
    pair_SNP_dict_with_location_all_genes = {}
    
    # major dict for GWASs
    # key:   GWAS_eQTL
    # value: SNPlist in the format of (GSNP_name,GSNP_chrom,GSNP_abs,GSNP_pval)
    GWAS_SNPlist_dict = {}
    # major dict for eQTLs
    #key:   GWAS_eQTL
    #value: gene_SNPlist_dict:  
            #key:   gene_name
            #value: SNPlist in the format of (eSNP_name,eSNP_chrom,eSNP_abs,eSNP_pval,aligned,tagged,gene,aligned_GSNP_name)                   
    eQTL_gene_SNPlist_dict = {}
   
    location_gene_dict = dao.Manhattan_gen_all_location_gene_dict()
 
    start_time = time.time()
    for gene in genes:

        # fetch pair_SNP_dict with only name and pval, then add location information
        #pair_SNP_dict = Manhattan_gen_pair_SNP_dict(web_disease_list,web_eQTL_list,pairNames,gene) 
        #pair_SNP_dict_with_location_curr_gene,GWAS_SNPlist_dict_curr_gene,eQTL_SNPlist_dict_curr_gene = dao.Manhattan_enhance_pair_SNP_dict_with_location(pair_SNP_dict,gene)
        
        # directly generate GWAS_SNPlist and eQTL_SNPlist for current gene across all GWAS_eQTL pairs 
        GWAS_SNPlist_dict_curr_gene,eQTL_SNPlist_dict_curr_gene = dao.fetch_pair_SNP_raw(web_disease_list,web_eQTL_list,gene,GSNP_cutoff)

        GWAS_SNPlist_dict = Manhattan_add_GWAS_SNPlist_dict(GWAS_SNPlist_dict,GWAS_SNPlist_dict_curr_gene)
        eQTL_gene_SNPlist_dict = Manhattan_add_eQTL_gene_SNPlist_dict(eQTL_gene_SNPlist_dict,eQTL_SNPlist_dict_curr_gene,gene) 

        #pair_SNP_dict_with_location_all_genes = Manhattan_add_pair_SNP_dict(pair_SNP_dict_with_location_all_genes,pair_SNP_dict_with_location_curr_gene)
    print 'fetching Manhattan SNP_list takes {} seconds'.format(time.time() - start_time)
    start_time = time.time()
    
    chrom_starts = dao.Manhattan_gen_chrom_starts()    

    eQTL_SNPlist_dict = dao.Manhattan_gen_eQTL_SNPlist(location_pval_chrom_SNPlist_dict,genes)

    gene_location_dict = dao.Manhattan_gen_gene_location_dict(genes) 

    all_eQTL_names = dao.Manhattan_get_all_eQTL_names()

    #location_pval_chrom_SNPlist_dict_1215, eQTL_SNPlist_dict_1215 = Manhattan_debug(pair_SNP_dict_with_location_all_genes)

    #Manhattan_pairNames = location_pval_chrom_SNPlist_dict_1215.keys()
    Manhattan_pairNames = GWAS_SNPlist_dict.keys()

    ret = {}
    #ret['location_pval_chrom_SNPlist_dict'] = location_pval_chrom_SNPlist_dict_1215
    ret['GWAS_SNPlist_dict'] = GWAS_SNPlist_dict
    ret['eQTL_gene_SNPlist_dict'] = eQTL_gene_SNPlist_dict
    ret['Manhattan_pairNames'] = list(Manhattan_pairNames)
    ret['Manhattan_geneNames'] = genes
    ret['chrom_starts'] = chrom_starts
    ret['gene_location_dict'] = gene_location_dict
    #ret['eQTL_SNPlist_dict'] = eQTL_SNPlist_dict_1215
    print 'Manhattan backend takes ' + str(time.time() - start_time_Manhattan) + ' seconds'
    return jsonify(ret)

# display only those GWAS SNPs that more significant than cutoff, for example 10-3
def Manhattan_filter_pair_SNP_dict_by_GWAS_pval_cutoff(pair_SNP_dict,Manhattan_GWAS_cutoff):
    pair_SNP_dict_filtered = {}
    for pair in pair_SNP_dict:
        orig_SNP_tuple_list = pair_SNP_dict[pair]
        new_SNP_tuple_list = []
        for SNP_tuple in orig_SNP_tuple_list:
            curr_GWAS_pval_str = SNP_tuple[2]   #GSNP_name,eSNP_name,Gpval,epval 
            curr_GWAS_pval = -1
            try:
                curr_GWAS_pval = float(curr_GWAS_pval_str)
            except ValueError:
                a = 1   #placeholder
            
            if curr_GWAS_pval < Manhattan_GWAS_cutoff:  # including those dummyone from pair_SNP_dict
                new_SNP_tuple_list.append(SNP_tuple)
            else:
                new_SNP_tuple_list.append(('dummy', 'dummy', '-1', '-1'))
        pair_SNP_dict_filtered[pair] = new_SNP_tuple_list
    return pair_SNP_dict_filtered 

def fetch_and_build_matrix(consider_all_genes_in_database):
    web_disease_list = session['web_disease_list']
    web_eQTL_list = session['web_eQTL_list']
    web_num_genes_per_pair = session['web_num_genes_per_pair']
    dao = getattr(g, 'dao', None)
    start_time = time.time()
    gene_p_qs,filtered_gene_names,gene_descriptions = dao.fetch_pair_gene(web_disease_list,web_eQTL_list,web_num_genes_per_pair,consider_all_genes_in_database)
    print "fetch_pair_gene take {} seconds".format(time.time() - start_time)
    if gene_p_qs is None:
        return None,None,None 
     
    return gene_p_qs,filtered_gene_names,gene_descriptions

@app.route('/')
def show_matrix():
    disease_names = sorted(disease_GWAS_dict.keys())
     
    if eQTL_name_all not in eQTL_names:
        eQTL_names.append(eQTL_name_all)
 
    if 'web_disease_list' not in session or 'web_eQTL_list' not in session:
        return render_template('show_matrix.html',eQTL_names = eQTL_names,disease_names = disease_names)
    
    try:
        page = int(request.args.get('page', 1))
    except ValueError:
        page = 1
    session['page'] = page

    gene_p_qs,filtered_gene_names,gene_descriptions = fetch_and_build_matrix(False) 
    if gene_p_qs is None:
        return render_template('show_matrix.html',eQTL_names = eQTL_names,disease_names = disease_names) 
    
    ret = {}
    ret['filtered_gene_names'] = filtered_gene_names
    ret['gene_descriptions'] = gene_descriptions
    ret['gene_p_qs'] = gene_p_qs
    ret['sorted_pair_names'] = sorted(gene_p_qs.keys())
    draw_pair_json_obj = json.dumps(ret)

    show_discover_sub_clusters_button = True
    #show_discover_sub_clusters_button = False
    #if len(gene_p_qs) <= 15:
    #    show_discover_sub_clusters_button = True   

    output_matrix_to_txt(ret)

    return render_template('show_matrix.html', pagination=None, page=page, eQTL_names=eQTL_names, disease_names = disease_names, draw_pair_json_obj=draw_pair_json_obj,show_discover_sub_clusters_button = show_discover_sub_clusters_button)

@app.route('/draw', methods=['POST'])
def draw():
    web_eQTL_list = ''

    if len(request.form.getlist(eQTL_name_all)) > 0:
        web_eQTL_list = ' '.join(eQTL_names) 
    else:
        for eQTL_name in eQTL_names:
            eQTL_name_selected_list = request.form.getlist(eQTL_name)
            if len(eQTL_name_selected_list) > 0:
                web_eQTL_list = web_eQTL_list + eQTL_name + ' ' 
           
    web_disease_list = request.form['disease_list']
    #web_eQTL_list  = request.form['eQTL_list']
    web_num_genes_per_pair = request.form['num_genes_per_pair']

    session['web_disease_list'] = web_disease_list
    session['web_eQTL_list'] = web_eQTL_list
    session['web_num_genes_per_pair'] = web_num_genes_per_pair
    return redirect(url_for('show_matrix'))

if __name__ == '__main__':
    #ip_for_current_machine = socket.gethostbyname(socket.gethostname())
    #app.run(host=ip_for_current_machine,port=55555,threaded=True)
    #app.run(host='localhost',port=15213,threaded=True)
#    app.run(host='169.230.81.176',port=55555,threaded=True)

    app.run(host='0.0.0.0',port=55556,threaded=True)

