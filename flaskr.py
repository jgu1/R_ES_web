#all the imports

import pdb
import MySQLdb
from flask import Flask, request, session, g, redirect, url_for, \
     abort, render_template, flash, jsonify
from contextlib import closing
from time import sleep
from flask.ext.paginate import Pagination
import re
import socket
from datetime import datetime,timedelta
from db_classes import DAO
import json
import pickle
# configuration
DEBUG = True
SECRET_KEY = 'development key'
eQTL_names = ['Dixon_07','Duan_08','Liang_2012','Muther_12','Myers_07','Schadt_08','Wright_14_pruned_e6_2','Zeller_10','merged_pickle']
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

@app.route('/detail')
@app.route("/data/<string:gene>")
def detail():
    gene = request.args.get('gene', 'empty')
    web_GWAS_list = session['web_GWAS_list']
    web_eQTL_list = session['web_eQTL_list']
    dao = getattr(g, 'dao', None)
    pair_SNP_dict,all_SNPs_list = dao.fetch_pair_SNP(web_GWAS_list,web_eQTL_list,gene)
    ret = {}
    ret['gene'] = gene
    ret['pair_SNP_dict'] = pair_SNP_dict
    ret['all_SNPs_list'] = all_SNPs_list
    return jsonify(ret)

def fetch_and_build_matrix():
    web_GWAS_list = session['web_GWAS_list']
    web_eQTL_list = session['web_eQTL_list']
    dao = getattr(g, 'dao', None)
    gene_p_qs,filtered_gene_names = dao.fetch_pair_gene(web_GWAS_list,web_eQTL_list)
    if gene_p_qs is None:
        return None,None,None 
     
    return gene_p_qs,None,filtered_gene_names

@app.route('/')
def show_matrix():
   
    if 'web_GWAS_list' not in session or 'web_eQTL_list' not in session:
        return render_template('show_matrix.html',eQTL_names = eQTL_names)
    
    try:
        page = int(request.args.get('page', 1))
    except ValueError:
        page = 1
    session['page'] = page

    gene_p_qs_for_this_page,pagination,filtered_gene_names_for_this_page = fetch_and_build_matrix() 
    if gene_p_qs_for_this_page is None:
        return render_template('show_matrix.html',eQTL_names = eQTL_names) 
    ret = {}
    ret['filtered_gene_names_for_this_page'] = filtered_gene_names_for_this_page
    ret['gene_p_qs_for_this_page'] = gene_p_qs_for_this_page
    ret['sorted_pair_names'] = sorted(gene_p_qs_for_this_page.keys())
    draw_pair_json_obj = json.dumps(ret)
    return render_template('show_matrix.html', pagination=pagination, page=page, eQTL_names=eQTL_names, draw_pair_json_obj=draw_pair_json_obj)

@app.route('/draw', methods=['POST'])
def draw():
   
    web_eQTL_list = ''
    for eQTL_name in eQTL_names:
        eQTL_name_selected_list = request.form.getlist(eQTL_name)
        if len(eQTL_name_selected_list) > 0:
            web_eQTL_list = web_eQTL_list + eQTL_name + ' ' 
    
    if not session.get('logged_in'):
        abort(401)
    
    web_GWAS_list = request.form['GWAS_list']
    #web_eQTL_list  = request.form['eQTL_list']

    session['web_GWAS_list'] = web_GWAS_list
    session['web_eQTL_list'] = web_eQTL_list
    return redirect(url_for('show_matrix'))

if __name__ == '__main__':
    ip_for_current_machine = socket.gethostbyname(socket.gethostname())
    #app.run(host=ip_for_current_machine,port=55555,threaded=True)
#    app.run(host='localhost',port=55555,threaded=True)
    app.run(host='169.230.81.176',port=55555,threaded=True)


