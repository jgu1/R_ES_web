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
@app.route("/sub_clusters/<string:alg><string:abs_cutoff><string:per_cutoff><string:converge_epsilon><string:converge_depth><string:est_col_width><string:filter_ratio><string:consider_all_genes_in_database>")
@app.route("/sub_clusters/<string:alg><string:consider_all_genes_in_database>")
def sub_clusters():

    alg = request.args.get('alg','')
    if alg == 'ISA':
        
        abs_cutoff       = request.args.get('abs_cutoff', '')
        per_cutoff       = request.args.get('per_cutoff','')
        converge_epsilon = request.args.get('converge_epsilon','')
        converge_depth   = request.args.get('converge_depth','')
        est_col_width    = request.args.get('est_col_width','')
        filter_ratio     = request.args.get('filter_ratio','')
        consider_all_genes_in_database = request.args.get('consider_all_genes_in_database','')
        if abs_cutoff == '':
            abs_cutoff = 3
        if per_cutoff == '':
            per_cutoff  = 0.5
        if converge_epsilon == '':
            converge_epsilon = 0.1
        if converge_depth == '':
            converge_depth  = 100
        if est_col_width == '':
            est_col_width = 20
        if filter_ratio == '':
            filter_ratio = 0.3
        if consider_all_genes_in_database == 'true':
            consider_all_genes_in_database = True
        else:
            consider_all_genes_in_database = False

        gene_p_qs,filtered_gene_names,gene_descriptions = fetch_and_build_matrix(consider_all_genes_in_database)
        sub_clusters = R_discover_sub_clusters_ISA(gene_p_qs,float(abs_cutoff),float(per_cutoff),float(converge_epsilon),float(converge_depth),float(est_col_width),float(filter_ratio))
   
    
    elif alg == 'PLAID':
       
        pdb.set_trace() 
        consider_all_genes_in_database = request.args.get('consider_all_genes_in_database','')
        if consider_all_genes_in_database == 'true':
            consider_all_genes_in_database = True
        else:
            consider_all_genes_in_database = False
        gene_p_qs,filtered_gene_names,gene_descriptions = fetch_and_build_matrix(consider_all_genes_in_database)
        sub_clusters = R_discover_sub_clusters_PLAID(gene_p_qs)
 
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

@app.route('/Manhattan')
@app.route("/Manhattan/<string:geneNames><string:pairNames>")
def Manhattan():
    geneNames = request.args.get('geneNames', 'empty')
    pairNames = request.args.get('pairNames','empty')

    web_disease_list = session['web_disease_list']
    web_eQTL_list = session['web_eQTL_list']

    genes = geneNames.split(',')
    genes = filter(None, genes) #remove empty gene, empty gene is unchecked in the web_interface   
    genes.sort()
 
    dao = getattr(g, 'dao', None)
    location_pval_chrom_SNPlist_dict = {}
    Manhattan_pairNames = set()
    start_time = time.time()
    for gene in genes:
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

        location_pval_chrom_SNPlist_dict_gene,Manhattan_pairNames_gene = dao.Manhattan_build_Manhattan_SNP_fields_list_dict(pair_SNP_dict,gene)
        for pairName in Manhattan_pairNames_gene:
            if pairName not in Manhattan_pairNames: # if this is a new pair
                Manhattan_pairNames.add(pairName)
                location_pval_chrom_SNPlist_dict[pairName] = location_pval_chrom_SNPlist_dict_gene[pairName]
                continue
            else:
                existing_list = location_pval_chrom_SNPlist_dict[pairName]
                new_list = location_pval_chrom_SNPlist_dict_gene[pairName] + existing_list
                location_pval_chrom_SNPlist_dict[pairName] = new_list 
    print 'fetching Manhattan SNP_list takes {} seconds'.format(time.time() - start_time)
    start_time = time.time()
    location_pval_chrom_SNPlist_dict,chrom_starts = dao.Manhattan_gen_abs_location_chrom(location_pval_chrom_SNPlist_dict)
    print 'generate abs_location and chrome takes {} seconds'.format(time.time() - start_time)
   
    gene_location_dict = dao.Manhattan_gen_gene_location_dict(genes) 

    ret = {}
    ret['location_pval_chrom_SNPlist_dict'] = location_pval_chrom_SNPlist_dict
    ret['Manhattan_pairNames'] = list(Manhattan_pairNames)
    ret['chrom_starts'] = chrom_starts
    ret['Manhattan_geneNames'] = genes
    ret['gene_location_dict'] = gene_location_dict
    return jsonify(ret)


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

