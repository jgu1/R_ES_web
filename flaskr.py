#all the imports

import pdb
import sqlite3
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

#app = Blueprint('papers',__name__)

def connect_db():
#    db =  sqlite3.connect(app.config['DATABASE'])
    '''
    db = MySQLdb.connect(host="localhost", # your host, usually localhost
        user="root", # your username
        passwd="genome", # your password
        db="ES_OUTPUT") # name of the data base
    return db
    '''
    return DAO()


'''
def init_db():
    with closing(connect_db()) as db:
        with app.open_resource('create_pubmed_cache.sql', mode='r') as f:
            db.cursor().executescript(f.read())
        db.commit()
'''

@app.before_request
def before_request():
    g.dao = connect_db()

@app.teardown_request
def teardown_request(exception):
    dao = getattr(g, 'dao', None)
    if dao is not None:
        dao.db.close()

def fetch_db_for_search_terms(disease,genes_included,genes_excluded):
    papers=[]
    count_dict={}
    disease = '+'.join(disease.split())
    for gene in genes_included:
        search_term = 'AND+' + disease + '+' + gene
        if not not genes_excluded:
            search_term += '+NOT+'+ '+'.join(genes_excluded)

        papers_for_curr_search_term = fetch_db_for_one_search_term(search_term)
        count_dict[gene]=len(papers_for_curr_search_term)
        papers += papers_for_curr_search_term
    return papers,count_dict

def fetch_db_for_one_search_term(search_term):
    sql_fetch_all_papers_for_search_term = ('with paper_ids as' 
         ' (select paper_id from search_terms as S ,term_paper_relation as T' 
         ' where S.search_term="'+ search_term +'" and S.id = T.term_id)' 
         ' select title,link,authors_str,journal_title,publish_time_str,abstract,keywords_str from paper_ids, papers' 
         ' where paper_ids.paper_id = papers.id;')
    cur = g.db.execute(sql_fetch_all_papers_for_search_term)
    papers = [dict(title=row[0], link=row[1], authors_str=row[2],journal_title=row[3],publish_time_str=row[4],abstract=row[5],keywords_str=row[6],search_term=search_term) for row in cur.fetchall()]
    return papers

def delete_db_for_one_search_term(search_term):
    sql_delete_all_papers_for_search_term = ('with paper_ids as' 
         ' (select paper_id from search_terms as S ,term_paper_relation as T' 
         ' where S.search_term="'+ search_term +'" and S.id = T.term_id)' 
         ' delete from papers' 
         ' where papers.id in paper_ids;')
    cur = g.db.execute(sql_delete_all_papers_for_search_term)
    
    sql_delete_all_term_paper_relations_for_search_term = ('with term_ids as'
        ' (select term_id from term_paper_relation as T, search_terms as S'
        ' where T.term_id = S.id and S.search_term="' + search_term + '") '
        ' delete from term_paper_relation'
        ' where term_paper_relation.term_id in term_ids ; ')
    cur = g.db.execute(sql_delete_all_term_paper_relations_for_search_term)

    sql_delete_search_terms_for_search_term = ('delete from search_terms'
        ' where search_terms.search_term="' + search_term + '" ;'
        )
    cur = g.db.execute(sql_delete_search_terms_for_search_term)
    g.db.commit()
'''
@app.route('/data')
@app.route("/data/<int:page>")
def data():

    page = session['page']
    dao = getattr(g, 'dao', None)
    gene_p_qs = dao.fetch_all_gene_p_q()
    gene_p_qs_for_this_page = {}
    for key in gene_p_qs:
         curr_rows = gene_p_qs[key]
         gene_p_qs_for_this_page[key] = curr_rows[(page-1)*GENE_P_Q_PER_PAGE:page*GENE_P_Q_PER_PAGE]
    return jsonify(gene_p_qs_for_this_page)
'''
@app.route('/detail')
@app.route("/data/<string:GWAS>/<string:eQTL>/<string:gene>")
def detail():
    GWAS = request.args.get('GWAS', 'empty')
    eQTL = request.args.get('eQTL', 'empty')
    gene = request.args.get('gene', 'empty')
    dao = getattr(g, 'dao', None)
    list_detail = dao.fetch_detail(GWAS,eQTL,gene)
    ret = {}
    ret['GWAS'] = GWAS
    ret['eQTL'] = eQTL
    ret['gene'] = gene
    ret['SNP_list'] = list_detail
    return jsonify(ret)
@app.route('/')
def show_papers():

    ''' 
    if 'disease' not in session or 'genes_included' not in session:
        return render_template('show_papers.html')

 
    disease = session['disease']
    genes_included = session['genes_included']
    genes_excluded = session['genes_excluded']

    all_papers,count_dict = fetch_db_for_search_terms(disease,genes_included,genes_excluded)
    '''
    dao = getattr(g, 'dao', None)
    #begin pagination 
    try:
        page = int(request.args.get('page', 1))
    except ValueError:
        page = 1
    
    session['page'] = page
    GENE_P_Q_PER_PAGE= app.config['GENE_P_Q_PER_PAGE'] 
    
    web_GWAS_list = session['web_GWAS_list']
    web_eQTL_list = session['web_eQTL_list']

    if not web_GWAS_list or not web_eQTL_list:
        return render_template('show_papers.html')
    
    gene_p_qs = dao.fetch_pair_gene(web_GWAS_list,web_eQTL_list)
    gene_p_qs_for_this_page = {}
    orig_length = -1
    for pair_name in gene_p_qs:
        orig_length_result = gene_p_qs[pair_name]
        if orig_length <0:
            orig_length = len(orig_length_result)   # get the length for pagination
        gene_p_qs_for_this_page[pair_name] = orig_length_result[(page-1)*GENE_P_Q_PER_PAGE:page*GENE_P_Q_PER_PAGE]       
 
    pagination = Pagination(page=page, total=orig_length, per_page=GENE_P_Q_PER_PAGE, record_name='genes for pairs')
    
    #end pagination
    gene_p_qs_json_obj = json.dumps(gene_p_qs_for_this_page)
     
    return render_template('show_papers.html',pagination=pagination, page=page,gene_p_qs_json_obj = gene_p_qs_json_obj)
    '''
    gene_p_qs_for_this_page = gene_p_qs[(page-1)*GENE_P_Q_PER_PAGE:page*GENE_P_Q_PER_PAGE] 
    pagination = Pagination(page=page, total=len(gene_p_qs), per_page=GENE_P_Q_PER_PAGE, record_name='gene_p_qs')
    #end pagination
    
    session['page'] = page
    return render_template('show_papers.html',gene_p_qs=gene_p_qs_for_this_page,pagination=pagination,page=page)
    '''
    
    '''
    gene_p_qs = dao.fetch_all_gene_p_q()
    keys = gene_p_qs.keys()
    pagination = Pagination(page=page, total=len(gene_p_qs[keys[0]]), per_page=GENE_P_Q_PER_PAGE, record_name='gene_p_qs')
    return render_template('show_papers.html',pagination=pagination,page=page)
    '''
def pop_db(disease,genes_included,genes_excluded):
    
    disease='+'.join(disease.split())
 
    for gene in genes_included:
        search_term = 'AND+' + disease + '+' + gene
        if not not genes_excluded:
            search_term += '+NOT+'+ '+'.join(genes_excluded)
        #pdb.set_trace()
        sql_check_if_search_term_has_results_already = ('select last_update from search_terms'
             ' where search_term="' + search_term +'"')
        cur = g.db.execute(sql_check_if_search_term_has_results_already)
        result = cur.fetchall()
        if len(result) > 0:
            time_delta = datetime.now() - datetime.strptime(result[0][0], "%Y-%m-%d %H:%M:%S.%f")
            if time_delta.seconds > 60*60*4: # re-fetch if the record is older than 1 day
                delete_db_for_one_search_term(search_term)
                esearch_fetch_parse.Main(DATABASE,search_term,gene)
        else:
            esearch_fetch_parse.Main(DATABASE,search_term,gene)

def parse_web_search_term(web_search_term_disease,web_search_term_genes_included,web_search_term_genes_excluded):
    disease = web_search_term_disease.strip()
    genes_included = web_search_term_genes_included.strip().split()
    genes_excluded = web_search_term_genes_excluded.strip().split()
    return disease,genes_included,genes_excluded

@app.route('/draw', methods=['POST'])
def draw():
 
    if not session.get('logged_in'):
        abort(401)
    
    web_GWAS_list = request.form['GWAS_list']
    web_eQTL_list  = request.form['eQTL_list']

    session['web_GWAS_list'] = web_GWAS_list
    session['web_eQTL_list'] = web_eQTL_list
    '''
    dao = getattr(g, 'dao', None) 
    gene_p_qs = dao.fetch_pair_gene(web_GWAS_list,web_eQTL_list)

    gene_p_qs_json_obj = json.dumps(gene_p_qs)
     
    pdb.set_trace()
    '''
    return redirect(url_for('show_papers'))

@app.route('/choose_term', methods=['GET'])
def choose_term():
    gene=request.args.get('gene', 1)
    session['genes_included']=[gene]
    return redirect(url_for('show_papers'))

@app.route('/login', methods=['GET', 'POST'])
def login():
    error = None
    if request.method == 'POST':
        if request.form['username'] not in app.config['USERNAME_PASSWORD_DICT']:
            error = 'Invalid username'
        else:
            session['logged_in'] = True
            flash('You were logged in')
            return redirect(url_for('show_papers'))
    return render_template('login.html', error=error)

@app.route('/logout')
def logout():
    session.pop('logged_in', None)
    flash('You were logged out')
    return redirect(url_for('show_papers'))

if __name__ == '__main__':
    ip_for_current_machine = socket.gethostbyname(socket.gethostname())
    #app.run(host=ip_for_current_machine,port=55555,threaded=True)
    app.run(host='localhost',port=55555,threaded=True)


