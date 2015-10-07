import os

abspath = os.path.abspath(__file__)
dirname = os.path.dirname(abspath)

csx_server_path = os.path.abspath(
    os.path.join(os.path.dirname( __file__ ), '..', 'csx_server')
)

csx_server_path += '/'

def writeHeaderHTML(path, version):
    html = path + "result_sheet.html"
    html = open(html, 'w')

    html.write("""
<!DOCTYPE html>
<html class="" lang="en">
<head>
  <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
  <title>CoNSENsX result sheet</title>
  <meta name="description" content="COmpliance of NMR Structural
                                    ENSembles with eXperimental data" />
  <link href="{0}public/_include/css/bootstrap.min.css" rel="stylesheet">
  <link href="{0}public/_include/css/style.css" rel="stylesheet">
  <link href="{0}public/_include/css/main.css" rel="stylesheet">
  <link href='http://fonts.googleapis.com/css?family=Titillium+Web:400,200,200italic,300,300italic,400italic,600,600italic,700,700italic,900' rel='stylesheet' type='text/css'>
</head>""".format(csx_server_path))

    html.close()


def writeFileTable(path, args, my_PDB, my_id, PDP_model_num,
                   my_NOE="not present", NOE_restraint_count=0):
    html = path + "result_sheet.html"
    html = open(html, 'a')

    html.write("""
<body>
<div id="shortcodes" class="page" style="padding: 30px 0 0 0;">
    <div class="container" style="width: 1240px;">
    <div class="row" id="res-page-row" style="margin-left: 0px;">
      <div class="span12">
        <div class="title-page">
        <h2 class="title">CoNSENsX</h2>
        <h3 class="title-description"><b class="red">Co</b>mpliance of <b class="red">
        N</b>MR-derived <b class="red">S</b>tructural <b class="red">Ens</b>embles with e<b class="red">x</b>perimental data</h3><br>
        <h3 class="title-description">Results sheet ID: <b class="red">{0}</b></h3>
        </div>
      </div>
    </div>

    <table class="files_table">
      <tr>
        <td class="head-td">PDB file:</td>
        <td><i>{1}</i></td>
        <td>{2} models found</td>
      </tr>
      <tr>
        <td class="head-td">NOE restraint file:</td>
        <td><i>{3}</i></td>
        <td>{4} distance restraints found</td>
      </tr>
      <tr>
        <td class="head-td">BMRB file:</td>
        <td><i>{5}</i></td>
        <td></td>
      </tr>
    </table>
    <div class="container">

    <button class="btn btn-primary" id="selection" type="button" data-toggle="collapse" data-target="#collapseExample" aria-expanded="false" aria-controls="collapseExample">
    Toggle selection options
    </button>
    <div class="collapse" id="collapseExample">
      <div>
        <h2 class="selection-heading">Compliance Measure</h2>
        <div class="btn-group group1" data-toggle="buttons">

          <label class="btn btn-primary">
            <input type="radio" name="measures" class="options" value="correlation" > correlation
          </label>

          <label class="btn btn-primary">
            <input type="radio" name="measures" class="options" value="q-value" > Q-value
          </label>

          <label class="btn btn-primary">
            <input type="radio" name="measures" class="options" value="rmsd" > RMSD
          </label>

        </div>
        <button class="btn btn-primary" id="startsel" type="button" onclick="myFunction()">
          Start selection
        </button>
      </div>
    </div>
  </div>""".format(my_id, my_PDB, PDP_model_num,
                       my_NOE, NOE_restraint_count, args.STR_file))

    html.close()


def writeRDC_table_open(path, name, RDC_list_num):
    html = path + "result_sheet.html"
    html = open(html, 'a')

    html.write("""
    <div class="results">
      <h4 class="table-source">{0} {1}</h4>""".format(name, RDC_list_num))
    html.close()


def writeRDC_table_close(path):
    html = path + "result_sheet.html"
    html = open(html, 'a')

    html.write("    </div>\n")
    html.close()


def writeRDC_data(path, RDC_type, used_values, correl, q_value, rmsd,
                  corr_graph_name, graph_name, mod_corr_graph_name,
                  input_id):
    html = path + "result_sheet.html"
    html = open(html, 'a')

    html.write("""
      <table class="result_table">
        <tbody>
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
          <td><img width="270" src="{7}"></td>
          </tr>
        </tbody>
      </table>\n
      <div class="involve involve-result">
      <form oninput="x.value=parseInt({8}.value)">
        <table class="selection-table">
          <tr>
            <td class="involve-table"><label>Involve with weight:</label></td>
            <td class="involve-table" style="width: 48px;"><output name="x" for="{8}">0</output></td>
            <td class="involve-table"><input type="range" class="inputrange" id="{8}" value="0" min="0" max="10"></td>
          </tr>
        </table>
      </form>
      </div>\n""".format(RDC_type, used_values,
                            '{0:.3f}'.format(correl),
                            '{0:.3f}'.format(q_value),
                            '{0:.3f}'.format(rmsd),
                             corr_graph_name, graph_name, mod_corr_graph_name,
                             input_id))

    html.close()


def write_table_open(path, table_name):
    html = path + "result_sheet.html"
    html = open(html, 'a')

    html.write("""
    <div class="results">
      <h4 class="table-source">{0}</h4>""".format(table_name))
    html.close()


def write_table_close(path):
    html = path + "result_sheet.html"
    html = open(html, 'a')

    html.write("    </div>\n")
    html.close()


def write_table_data(path, data_type, used_values, correl, q_value, rmsd,
                     corr_graph_name, graph_name, input_id,
                     mod_corr_graph_name=None):
    html = path + "result_sheet.html"
    html = open(html, 'a')

    if mod_corr_graph_name:
      html.write("""
        <table class="result_table">
          <tbody>
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
            <td><img width="270" src="{7}"></td>
            </tr>
          </tbody>
        </table>\n
        <div class="involve involve-result">
        <form oninput="x.value=parseInt({8}.value)">
          <table class="selection-table">
            <tr>
              <td class="involve-table"><label>Involve with weight:</label></td>
              <td class="involve-table" style="width: 48px;"><output name="x" for="{8}">0</output></td>
              <td class="involve-table"><input type="range" class="inputrange" id="{8}" value="0" min="0" max="10"></td>
            </tr>
          </table>
        </form>
        </div>\n""".format(data_type, used_values,
                              '{0:.3f}'.format(correl),
                              '{0:.3f}'.format(q_value),
                              '{0:.3f}'.format(rmsd),
                               corr_graph_name, graph_name,
                               mod_corr_graph_name, input_id))
    else:
        html.write("""
        <table class="result_table">
          <tbody>
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
            </tr>
          </tbody>
        </table>\n
        <div class="involve involve-result">
        <form oninput="x.value=parseInt({7}.value)">
          <table class="selection-table">
            <tr>
              <td class="involve-table"><label>Involve with weight:</label></td>
              <td class="involve-table" style="width: 48px;"><output name="x" for="{7}">0</output></td>
              <td class="involve-table"><input type="range" class="inputrange" id="{7}" value="0" min="0" max="10"></td>
            </tr>
          </table>
        </form>
        </div>\n""".format(data_type, used_values,
                              '{0:.3f}'.format(correl),
                              '{0:.3f}'.format(q_value),
                              '{0:.3f}'.format(rmsd),
                               corr_graph_name, graph_name, input_id))

    html.close()


def write_bottom_table(path, NOE_violations, PRIDE_data):
    html = path + "result_sheet.html"
    html = open(html, 'a')

    html.write("""
    <div class="results">
      <h4 class="table-source">NOE violations and PRIDE-NMR</h4>
      <table class="result_table">

          <td><table class="NOE_PRIME_table">
          <tr><td><strong>NOE distance violation</strong></td></tr>
          <tr><td>Total # of violations:</td><td>{0}</td></tr>
          </table></td>
          <td width="10"></td>
          <td><img width="320" src="NOE_hist.svg"></td>

          <td width="30"></td>

          <td><table class="NOE_PRIME_table">
          <tr><td><strong>PRIDE-NMR</strong></td></tr>
          <tr><td>Model with best score:</td><td>{1}</td></tr>
          <tr><td>Model with worst score:</td><td>{2}</td></tr>
          <tr><td>Average score:</td><td>{3}</td></tr>
          <tr><td>Standard Deviation:</td><td>{4}</td></tr>
          </table></td>
          <td width="10"></td>
          <td><img width="320" src="PRIDE-NMR_score.svg"></td>

      </table>
    </div>\n""".format(NOE_violations,
                       PRIDE_data[0],
                       PRIDE_data[1],
                       '{0:.3f}'.format(PRIDE_data[2]),
                       '{0:.3f}'.format(PRIDE_data[3])))

    html.close()


def close_HTML(path):
    html = path + "result_sheet.html"
    html = open(html, 'a')

    html.write("""
  </div>
</div>
<script src="{0}public/_include/js/jquery.js"></script>
<script src="{0}public/_include/js/bootstrap.min.js"></script>
<script src="{0}public/_include/js/involve.js"></script>
</body>
</html>""".format(csx_server_path))

    html.close()
