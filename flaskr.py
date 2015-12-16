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
# configuration
DEBUG = True
SECRET_KEY = 'development key'
USERNAME_PASSWORD_DICT={'hao':'genome','jiashun':'genome','erxin':'genome','jun':'genome','yanqiu':'genome','jialiang':'genome'}
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

@app.route('/sortAlongPair')
@app.route("/sortAlongPair/<string:sortAlongPairName>")
def sortAlongPair():
    sortAlongPairName = request.args.get('sortAlongPairName','empty')
    web_GWAS_list = session['web_GWAS_list']
    web_eQTL_list = session['web_eQTL_list']
    dao = getattr(g, 'dao', None)
    gene_p_qs,filtered_gene_names = dao.fetch_pair_gene(web_GWAS_list,web_eQTL_list)
    
    sortAlongPairList = gene_p_qs[sortAlongPairName]
    #change every dummy value to 1
    sortAlongPairListPvals = []
    for i in range(len(sortAlongPairList)):
        curr_gene_tuple = sortAlongPairList[i]
        if 'dummy' in curr_gene_tuple[0]:
            sortAlongPairListPvals.append(1.1)    # if empty, append a large value so it's moved back, real pval never exceed 1  
        else:
            sortAlongPairListPvals.append(float(curr_gene_tuple[3]))

    sort_idx = sorted(range(len(sortAlongPairListPvals)),key=lambda x:sortAlongPairListPvals[x])
    filtered_gene_names = [filtered_gene_names[i] for i in sort_idx]
     
    for pair_name in gene_p_qs:
        orig_order_list = gene_p_qs[pair_name]
        sort_order_list = [orig_order_list[i] for i in sort_idx]   
        gene_p_qs[pair_name] = sort_order_list
 
    gene_p_qs_for_this_page = {}
    page = 1
    max_length = -1
    for pair_name in gene_p_qs:
        orig_length_result = gene_p_qs[pair_name]
        if max_length < len(orig_length_result):
            max_length = len(orig_length_result)   # get the length for pagination
        gene_p_qs_for_this_page[pair_name] = orig_length_result[(page-1)*GENE_P_Q_PER_PAGE:page*GENE_P_Q_PER_PAGE]       
 
    pagination = Pagination(page=page, total=max_length, per_page=GENE_P_Q_PER_PAGE, record_name='genes for pairs')
 
    #end pagination
    filtered_gene_names_for_this_page = filtered_gene_names[(page-1)*GENE_P_Q_PER_PAGE:page*GENE_P_Q_PER_PAGE]
    ret = {}
    ret['filtered_gene_names_for_this_page'] = filtered_gene_names_for_this_page
    ret['gene_p_qs_for_this_page'] = gene_p_qs_for_this_page
    ret['sorted_pair_names'] = sorted(gene_p_qs_for_this_page.keys())
    return jsonify(ret)  

@app.route('/')
def show_matrix():
    dao = getattr(g, 'dao', None)
    #begin pagination 
    try:
        page = int(request.args.get('page', 1))
    except ValueError:
        page = 1
    
    session['page'] = page
    GENE_P_Q_PER_PAGE= app.config['GENE_P_Q_PER_PAGE'] 
    if 'web_GWAS_list' not in session or 'web_eQTL_list' not in session:
        return render_template('show_matrix.html')

    web_GWAS_list = session['web_GWAS_list']
    web_eQTL_list = session['web_eQTL_list']

    gene_p_qs,filtered_gene_names = dao.fetch_pair_gene(web_GWAS_list,web_eQTL_list)
    gene_p_qs_for_this_page = {}
    max_length = -1
    for pair_name in gene_p_qs:
        orig_length_result = gene_p_qs[pair_name]
        if max_length < len(orig_length_result):
            max_length = len(orig_length_result)   # get the max length for pagination
        gene_p_qs_for_this_page[pair_name] = orig_length_result[(page-1)*GENE_P_Q_PER_PAGE:page*GENE_P_Q_PER_PAGE]       
         
    pagination = Pagination(page=page, total=max_length, per_page=GENE_P_Q_PER_PAGE, record_name='genes for pairs')
   
    #end pagination
    filtered_gene_names_for_this_page = filtered_gene_names[(page-1)*GENE_P_Q_PER_PAGE:page*GENE_P_Q_PER_PAGE]
    ret = {}
    ret['filtered_gene_names_for_this_page'] = filtered_gene_names_for_this_page
    ret['gene_p_qs_for_this_page'] = gene_p_qs_for_this_page
    ret['sorted_pair_names'] = sorted(gene_p_qs_for_this_page.keys())
    draw_pair_json_obj = json.dumps(ret)
    return render_template('show_matrix.html',pagination=pagination, page=page,draw_pair_json_obj=draw_pair_json_obj)

@app.route('/draw', methods=['POST'])
def draw():
 
    if not session.get('logged_in'):
        abort(401)
    
    web_GWAS_list = request.form['GWAS_list']
    web_eQTL_list  = request.form['eQTL_list']

    session['web_GWAS_list'] = web_GWAS_list
    session['web_eQTL_list'] = web_eQTL_list
    return redirect(url_for('show_matrix'))

@app.route('/login', methods=['GET', 'POST'])
def login():
    error = None
    if request.method == 'POST':
        if request.form['username'] not in app.config['USERNAME_PASSWORD_DICT']:
            error = 'Invalid username'
        else:
            session['logged_in'] = True
            flash('You were logged in')
            return redirect(url_for('show_matrix'))
    return render_template('login.html', error=error)

@app.route('/logout')
def logout():
    session.pop('logged_in', None)
    flash('You were logged out')
    return redirect(url_for('show_matrix'))

if __name__ == '__main__':
    ip_for_current_machine = socket.gethostbyname(socket.gethostname())
    #app.run(host=ip_for_current_machine,port=55555,threaded=True)
    app.run(host='localhost',port=55555,threaded=True)


