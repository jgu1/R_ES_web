#all the imports

import pdb
import MySQLdb
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
from inspect_matrix import discover_sub_clusters
# configuration
DEBUG = True
SECRET_KEY = 'development key'
eQTL_names = ['Dixon_07','Duan_08','Liang_2012','Muther_12','Myers_07','Schadt_08','Wright_14_pruned_e6_2','Zeller_10','merged_pickle']
disease_GWAS_dict = pickle.load(open("disease_GWAS_dict.pickle","rb"))
GENE_P_Q_PER_PAGE=30
page=1

app = Flask(__name__)
app.config.from_object(__name__)


def connect_db():
    return DAO()

@app.before_request
def before_request():
    g.dao = connect_db()

@app.teardown_request
def teardown_request(exception):
    dao = getattr(g, 'dao', None)
    if dao is not None:
        dao.db.close()

@app.route('/sub_clusters')
def sub_clusters():
    gene_p_qs,pagination,filtered_gene_names,gene_descriptions = fetch_and_build_matrix()
    sub_clusters = discover_sub_clusters(gene_p_qs)
    
    serisables = []
    for sub_cluster in sub_clusters:
        row_comb = sub_cluster.row_comb
        cols = list(sub_cluster.cols)
        cols.sort()
        curr_serisable = [row_comb,cols]
        serisables.append(curr_serisable)

    ret = {}
    ret['gene_p_qs'] = gene_p_qs
    ret['gene_descriptions'] = gene_descriptions
    ret['pairname_idx'] = serisables

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
            pair_SNP_dict[pair] = pair_SNP_dict_all[pair]
    else:
        pair_SNP_dict = pair_SNP_dict_all
    ret = {}
    ret['gene'] = gene
    ret['pair_SNP_dict'] = pair_SNP_dict
    ret['all_SNPs_list'] = all_SNPs_list
    return jsonify(ret)

def fetch_and_build_matrix():
    web_disease_list = session['web_disease_list']
    web_eQTL_list = session['web_eQTL_list']
    dao = getattr(g, 'dao', None)
    gene_p_qs,filtered_gene_names,gene_descriptions = dao.fetch_pair_gene(web_disease_list,web_eQTL_list)
    if gene_p_qs is None:
        return None,None,None 
     
    return gene_p_qs,None,filtered_gene_names,gene_descriptions

@app.route('/')
def show_matrix():
    disease_names = sorted(disease_GWAS_dict.keys())
     
 
    if 'web_disease_list' not in session or 'web_eQTL_list' not in session:
        return render_template('show_matrix.html',eQTL_names = eQTL_names,disease_names = disease_names)
    
    try:
        page = int(request.args.get('page', 1))
    except ValueError:
        page = 1
    session['page'] = page

    gene_p_qs,pagination,filtered_gene_names,gene_descriptions = fetch_and_build_matrix() 
    if gene_p_qs is None:
        return render_template('show_matrix.html',eQTL_names = eQTL_names,disease_names = disease_names) 
    
    ret = {}
    ret['filtered_gene_names'] = filtered_gene_names
    ret['gene_descriptions'] = gene_descriptions
    ret['gene_p_qs'] = gene_p_qs
    ret['sorted_pair_names'] = sorted(gene_p_qs.keys())
    draw_pair_json_obj = json.dumps(ret)
    return render_template('show_matrix.html', pagination=pagination, page=page, eQTL_names=eQTL_names, disease_names = disease_names, draw_pair_json_obj=draw_pair_json_obj)

@app.route('/draw', methods=['POST'])
def draw():
    web_eQTL_list = ''
    for eQTL_name in eQTL_names:
        eQTL_name_selected_list = request.form.getlist(eQTL_name)
        if len(eQTL_name_selected_list) > 0:
            web_eQTL_list = web_eQTL_list + eQTL_name + ' ' 
       
    web_disease_list = request.form['disease_list']
    #web_eQTL_list  = request.form['eQTL_list']

    session['web_disease_list'] = web_disease_list
    session['web_eQTL_list'] = web_eQTL_list
    return redirect(url_for('show_matrix'))

if __name__ == '__main__':
    ip_for_current_machine = socket.gethostbyname(socket.gethostname())
    #app.run(host=ip_for_current_machine,port=55555,threaded=True)
#    app.run(host='localhost',port=55555,threaded=True)
    app.run(host='169.230.81.176',port=55555,threaded=True)


