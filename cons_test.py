#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import time
import os

import nmrpystar

version = "1.0"
pales   = "/home/daniel/Dokumente/önlab/gz_pack/pales/linux/pales"

#--------------------  Setting up and parsing arguments   --------------------#

consensx_usage = '%(prog)s -b STR_file -f PDB_file [options]'

parser = argparse.ArgumentParser(add_help = False,
                                 usage    = consensx_usage)

parser.add_argument("-b", "--STR_file", help = "restraints file")
parser.add_argument("-f", "--PDB_file", help = "PDB file (fitted)")
parser.add_argument("-r", "--XPLOR_file", default="",
                    help = "X-PLOR restraint file")

parser.add_argument("-h", "--help", help="show help and exit",
                    action='store_true')
parser.add_argument('-c', action='store_true',
                    help="show credits")

parser.add_argument("-l", "--lc_model", choices=["bic", "pf1"], default="bic",
                    help="rdc lc model: <bic|pf1>")
parser.add_argument("-R", action='store_true', default=False,
                    help="causes to do SVD for back-calculating RDC data")
parser.add_argument("-O", "--output_format", choices=["txt", "html"],
                    default="html|txt", help="output format")
parser.add_argument("-o", "--output_file_name", default="consensx",
                    help="output file name")
parser.add_argument("-T", "--line_thickness", help="line thickness in raphs")

args = parser.parse_args()


help_text = """
usage: cons_test.py -b <STR_file> -f <PDB_file>
             [-l <RDC_LC_MODEL: <bic|pf1>] [-R] [-O <txt|html>] [-o <outfile>]
             [-T <linethickness_in_graphs>]

  -h, --help            show this help message and exit
  -c                    show credits

  -b STR_file           restraints file
  -f PDB_file           PDB file (fitted)

  -l {bic,pf1}          rdc lc model: <bic|pf1>
  -R                    causes to do SVD for back-calculating RDC data
  -r                    X-PLOR restraint file
  -O {txt,html}         output format
  -o OUTPUT_FILE_NAME   output file name
  -T LINE_THICKNESS     line thickness in raphs


CoNSEnsX: assessing the compliance of varios NMR data with a protein
          structural ensemble
The program invokes SHIFTX for chemical shift and PALES for RDC calculation
See their respective documentation for more.
"""

cred = """
* Description of the CoNSEnsX concept and method can be found in:
-----------------------------------------------------------------
  CoNSEnsX: assessing the accuracy of NMR-derived protein structural ensembles
  Ángyán et al, 2009, submitted

* SHIFTX is described in:
------------------------
   Rapid and accurate calculation of protein 1H. 13C amd 15N
   Neal et al. (2003) J. Biomol. NMR 26:215.

* PALES is described in:
------------------------
  Prediction of sterically induced alignment in a dilute liquid
  crystalline phase: aid to protein strcuture determination by NMR
  Zweckstetter & Bax (2000) J. Am. Chem. Soc. 122:3791
"""

if args.c:
    print(cred)
    raise SystemExit

if args.help:
    print(help_text)
    raise SystemExit

if not args.STR_file or not args.PDB_file:
    print("missing file")
    parser.print_usage()
    raise SystemExit


#-------------------------  Setting up output files   ------------------------#
date = time.strftime("%a %d. %b %X %Z %Y")

if "txt" in args.output_format:
    txt = args.output_file_name + ".txt"
    txt = open(txt, 'w')

    txt.write("CoNSEnsX version " + version + " started on " + date + "\n")
    txt.write("========================================================\n")
    txt.write("Input files specified:" +
                   "\n\tPDB file: "                 + args.PDB_file +
                   "\n\tX-PLOR restraint file: "    + args.XPLOR_file +
                   "\n\tBMRB file: "                + args.STR_file + "\n\n")

if "html" in args.output_format:
    html = args.output_file_name + ".html"
    html  = open(html, 'w')

    html.write( "<html><head><title>CoNSEnsX output</title></head>\n" +
                "<body bgcolor=white><h2>CoNSEnsX version " + version +
                "</h2>\n<small>started on " + date + " </small><br>\n")

    html.write( "<br><b>Input files specified:</b>\n<ul>\n" +
                "<li>PDB file: " + args.PDB_file + "</li>\n" +
                "<li>X-PLOR restraint file: " + args.XPLOR_file + "</li>\n" +
                "<li>BMRB file: " + args.STR_file + "</li>\n</ul>\n<table>\n")


#-------------------------  Making temporary folder   ------------------------#
if not os.path.exists("temp"):
    os.makedirs("temp")

#--------------------------  Splitting of PDB file   -------------------------#

def pdb_splitter(PDB_file):
    try:
        my_pdb = open(PDB_file)
    except FileNotFoundError:
        print(PDB_file + " was not found")
        raise SystemExit

    model_names = []
    model_data  = []

    my_name = ""
    my_data = []

    for line in my_pdb:
        if line.startswith("MODEL"):
            my_name = line.strip().split()[1]
        elif line.startswith("ATOM") or line.startswith("TER"):
            my_data.append(line.strip())
        elif line.startswith("ENDMDL"):
            model_names.append(my_name)
            model_data.append(my_data)
            my_name = ""
            my_data = []
        else:
            continue

    for i in range(len(model_names)):
        file_name = "temp/model_" + model_names[i] + ".pdb"
        temp_pdb  = open(file_name, 'w')
        temp_pdb.write("HEADER    MODEL " + model_names[i] + "\n")

        for _ in model_data[i]:
            temp_pdb.write(_ + "\n")

        temp_pdb.write("END")
        temp_pdb.close()

pdb_splitter(args.PDB_file)


#------------------------------  Read STR file   -----------------------------#

star_file = open(args.STR_file)
myString = ""
for line in star_file:
    myString += line

#----------------------  Parse STR file with nmrpystar  ----------------------#

print "parsing " + args.STR_file + "..."
parsed = nmrpystar.parse(myString)
if parsed.status == 'success':
    print 'it worked!!  ' #, parsed.value
else:
    print 'uh-oh, there was a problem with the string I gave it ... ', parsed


def get_RDC_list(dataBlock, saveShiftName='RDC_list_1'):
    saveShifts = dataBlock.saves[saveShiftName]
    loopShifts = saveShifts.loops[-1]
    for ix in range(len(loopShifts.rows)):
        row = loopShifts.getRowAsDict(ix)
        # print row.keys()
        # print row.values(), "\n"

get_RDC_list(parsed.value)

for key in parsed.value.__dict__['saves'].keys():
    if key.startswith("RDC_list_"):
        print(key)
