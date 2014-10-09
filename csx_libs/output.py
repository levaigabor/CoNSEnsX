#!/usr/bin/python
# -*- coding: utf-8 -*-

import time
import os

date    = time.strftime("%a %d. %b %X %Z %Y")


def writeHeaderTXT(path, args, version):
    txt = path +"result_sheet.txt"
    txt = open(txt, 'w')

    txt.write("CoNSEnsX version " + version + " started on " + date + "\n")
    txt.write("========================================================\n")
    txt.write("Input files specified:" +
                   "\n\tPDB file: "                 + args.PDB_file +
                   "\n\tX-PLOR restraint file: "    + args.XPLOR_file +
                   "\n\tBMRB file: "                + args.STR_file + "\n\n")


def writeHeaderHTML(path, version):
    html = path + "result_sheet.html"
    html = open(html, 'w')

    html.write("""
<!DOCTYPE html>

<head>
  <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
  <title>CoNSENsX result sheet</title>
  <meta name="description" content="COmpliance of NMR Structural
                                    ENSembles with eXperimental data" />
  <link href="../../csx_ClientSide/_include/css/bootstrap.min.css" rel="stylesheet">
  <link href="../../csx_ClientSide/_include/css/main.css" rel="stylesheet">
  <link href="../../csx_ClientSide/_include/css/supersized.css" rel="stylesheet">
  <link href="../../csx_ClientSide/_include/css/bootstrap-responsive.min.css" rel="stylesheet">
  <link href="../../csx_ClientSide/_include/css/responsive.css" rel="stylesheet">
  <link href='http://fonts.googleapis.com/css?family=Titillium+Web:400,200,
              200italic,300,300italic,400italic,600,600italic,700,700italic,900'
              rel='stylesheet' type='text/css'>
</head>""")


def writeFileTable(path, args, my_PDB, my_id, PDP_model_num):
    html = path + "result_sheet.html"
    html = open(html, 'a')

    html.write("""
<body>
<div id="shortcodes" class="page" style="padding: 30px 0 0 0;">
    <div class="container">
    <div class="row" id="res-page-row">
      <div class="span12">
        <div class="title-page">
        <h2 class="title">CoNSENsX</h2>
        <h3 class="title-description"><b><font>Co</font></b>mpliance of <b>
        <font>N</font></b>MR-derived <b><font>S</font></b>tructural <b><font>Ens</font></b>embles with e<b><font>x</font></b>perimental data</h3><br>
        <h3 class="title-description">Results sheet ID: <font><b>{0}</b></font></h3>
        </div>
      </div>
    </div>

    <table class="files_table">
      <tr>
        <td id="head-td">PDB file:</td>
        <td><i>{1}</i></td>
        <td>{2} models found</td>
      </tr>
      <tr>
        <td id="head-td">X-PLOR restraint file:</td>
        <td><i></i></td>
        <td>0 distance restdaints found</td>
      </tr>
      <tr>
        <td id="head-td">BMRB file:</td>
        <td><i>{3}</i></td>
        <td></td>
      </tr>
    </table>""".format(my_id, my_PDB, PDP_model_num, args.STR_file) )


def writeRDC_table_open(path, name, RDC_list_num):
    html = path + "result_sheet.html"
    html = open(html, 'a')

    html.write("""
    <div class="results">
      <h4 class="table-source">{0} {1}</h4>
      <table class="result_table">\n""".format(name, RDC_list_num))


def writeRDC_table_close(path):
    html = path + "result_sheet.html"
    html = open(html, 'a')

    html.write("      </table>\n    </div>\n")


def writeRDC_data(path, RDC_type, used_values,
                  correl, q_value, rmsd, corr_graph_name, graph_name):

    html = path + "result_sheet.html"
    html = open(html, 'a')

    html.write("""
          <tr>
          <td><table class="values_table">
          <tr><td><strong>{0}</strong></td><td></td></tr>
          <tr><td>Values:</td><td>{1}</td></tr>
          <tr><td>Correlation:</td><td>{2}</td></tr>
          <tr><td>Q-factor:</td><td>{3} %</td></tr>
          <tr><td>RMSD:</td><td>{4}</td></tr>
          </table></td>
          <td><img width="270" src="{5}"></td>
          <td><img width="450" src="{6}"></td>
          </tr>\n""".format(RDC_type, used_values,
                            '{0:.3f}'.format(correl),
                            '{0:.3f}'.format(q_value),
                            '{0:.3f}'.format(rmsd),
                             corr_graph_name, graph_name))


def write_table_open(path, table_name):
    html = path + "result_sheet.html"
    html = open(html, 'a')

    html.write("""
    <div class="results">
      <h4 class="table-source">{0}</h4>
      <table class="result_table">\n""".format(table_name))


def write_table_close(path):
    html = path + "result_sheet.html"
    html = open(html, 'a')

    html.write("      </table>\n    </div>\n")


def write_table_data(path, data_type, used_values,
                     correl, q_value, rmsd, corr_graph_name, graph_name):
    html = path + "result_sheet.html"
    html = open(html, 'a')

    html.write("""
          <tr>
          <td><table class="values_table">
          <tr><td><strong>{0}</strong></td><td></td></tr>
          <tr><td>Values:</td><td>{1}</td></tr>
          <tr><td>Correlation:</td><td>{2}</td></tr>
          <tr><td>Q-factor:</td><td>{3} %</td></tr>
          <tr><td>RMSD:</td><td>{4}</td></tr>
          </table></td>
          <td><img width="270" src="{5}"></td>
          <td><img width="450" src="{6}"></td>
          </tr>\n""".format(data_type, used_values,
                            '{0:.3f}'.format(correl),
                            '{0:.3f}'.format(q_value),
                            '{0:.3f}'.format(rmsd),
                             corr_graph_name, graph_name))
