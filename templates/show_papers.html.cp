{% extends "layout.html" %}
{% block body %}
  <script src={{ url_for ('static',filename='jquery.min.js')}}></script>
  <script src={{ url_for ('static',filename='main.js')}}></script>
  <script src={{ url_for ('static',filename='d3.v3.min.js')}}></script>
  {% if session.logged_in %}
    <form action="{{ url_for('search') }}" method=post class=search>
      <dl>
        <dt>Disease/Trait:
        <dd><input type=text size=30 name=disease>
        <br>genes included<dd><input type=text size=30 name=genes_included>
        <br>genes excluded<dd><input type=text size=30 name=genes_excluded>
        <dd><input type=submit value=Search>
      </dl>
    </form>


    {{ pagination.info }}
    {{ pagination.links }} 
      {% for gene_p_q in gene_p_qs %}
      <ul class=gene_p_q>
        <li><span>GWAS: </span>{{gene_p_q.GWAS|safe}}<br>
        <li><span>eQTL: </span>{{gene_p_q.eQTL|safe}}<br>
        <li><span>gene: </span>{{gene_p_q.gene|safe}}<br>
        <li><span>pval: </span>{{gene_p_q.pval|safe}}<br>
        <li><span>qval: </span>{{gene_p_q.qval|safe}}<br>
      </ul>
      {% endfor%}
    {{ pagination.links }}
  {% endif %}
{% endblock %}