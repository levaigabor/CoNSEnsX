#!/usr/bin/python
# -*- coding: utf-8 -*-

# standard modules
from __future__ import print_function
import time
import sys
import os
import subprocess
import re
import math



# installed modules
import nmrpystar

# own modules
import csx_libs.misc    as csx_misc
import csx_libs.methods as csx_func
import csx_libs.objects as csx_obj

version = "1.0"


parser = csx_misc.createParser()
args = parser.parse_args()


if args.c:
    print(csx_misc.cred)
    raise SystemExit

if args.help:
    print(csx_misc.help_text)
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
if not os.path.exists("temp"):          # create temp folder
    os.makedirs("temp")
else:
    for f in os.listdir("temp"):        # clean temp folder
        os.remove("temp/" + f)

csx_func.pdb_cleaner(args.PDB_file)      # bringing PDB to format
csx_func.pdb_splitter(args.PDB_file)     # splitting of PDB file


#------------------------  Read  and parse STR file   ------------------------#

star_file = open(args.STR_file)
myString = ""

for line in star_file:
    myString += line

parsed = nmrpystar.parse(myString)

if parsed.status != 'success':
    print('Error during STR parsing: ', parsed)
    raise SystemExit


RDC_lists  = csx_func.get_RDC_lists(parsed.value)
pdb_models = os.listdir("temp")

##### -------------------------- for RDC_list in RDC_lists:

csx_func.callPalesOn(pdb_models, RDC_lists[0], args.lc_model, args.R)

averageRDC = csx_func.avgPalesRDCs("pales.out")


print("Correl: ", csx_func.calcCorrel(averageRDC, RDC_lists[0]))
print("Q-val:  ", csx_func.calcQValue(averageRDC, RDC_lists[0]))
print("RMSD:   ", csx_func.calcRMSD(averageRDC, RDC_lists[0]))
csx_func.makeGraph(averageRDC, RDC_lists[0])
csx_func.makeCorrelGraph(averageRDC, RDC_lists[0])

##### -------------------------- for RDC_list in RDC_lists


#-------------------------  S2 parse from STR file   -------------------------#
try:
    saveShifts    = parsed.value.saves["order_param"]
except KeyError:
    print("No S2 parameter list found")

loopShifts     = saveShifts.loops[-1]
S2_records = []

for ix in range(len(loopShifts.rows)):   # fetch values from STR file
            row = loopShifts.getRowAsDict(ix)

            S2_records.append(csx_obj.S2_Record(row["Residue_seq_code"],
                                                row["Atom_name"],
                                                row["S2_value"]))



#---------------------------------  S2 calc  ---------------------------------#


S2_calced = csx_func.calcS2(args.PDB_file, S2_records)

print("S2_corr:", csx_func.calcCorrel(S2_calced, S2_records))
print("S2Q-val:", csx_func.calcQValue(S2_calced, S2_records))
print("RMSD:   ", csx_func.calcRMSD(S2_calced, S2_records))
csx_func.makeGraph(S2_calced, S2_records)
csx_func.makeCorrelGraph(S2_calced, S2_records)

