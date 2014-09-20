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
import prody


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

#####for RDC_list in RDC_lists:
csx_func.callPalesOn(pdb_models, RDC_lists[0], args.lc_model, args.R)

averageRDC = csx_func.avgPalesRDCs("pales.out")

for RDC_type in averageRDC.keys():
    print(RDC_type)

print("Correl: ", csx_func.calcCorrel(averageRDC, '0_N_H', RDC_lists[0]))
print("Q-val:  ", csx_func.calcQValue(averageRDC, '0_N_H', RDC_lists[0]))
print("RMSD:   ", csx_func.calcRMSD(averageRDC, '0_N_H', RDC_lists[0]))
# csx_func.makeGraph(averageRDC, '0_N_H', RDC_lists[0])



#-------------------------  S2 parse from STR file   -------------------------#

class S2_Record(object):
    """Class for storing S2 data"""

    def __init__(self, resnum, S2_type, S2_value):
        self.resnum   = resnum
        self.S2_type  = S2_type
        self.S2_value = S2_value


try:
    saveShifts    = parsed.value.saves["order_param"]
except KeyError:
    print("No S2 parameter list found")

loopShifts     = saveShifts.loops[-1]
S2_records = []

for ix in range(len(loopShifts.rows)):   # fetch values from STR file
            row = loopShifts.getRowAsDict(ix)

            S2_records.append(S2_Record(row["Residue_seq_code"],
                                        row["Atom_name"],
                                        row["S2_value"]))




#-------------------------  S2 calc from PDB file   --------------------------#

# parsing PDB file into models (model_list)
model_list = []
model_num = 1

while True:
    try:
        model_list.append(prody.parsePDB(args.PDB_file, model=model_num, ter=True))
        model_num += 1
    except prody.proteins.pdbfile.PDBParseError:
        break


# get NH vectors from models (model_data[] -> vectors{resnum : vector})
model_data = []

for model in model_list:
    current_Resindex = 1
    has_H, has_N = False, False
    vectors = {}

    for atom in model:
        atom_res = atom.getResindex() + 1

        if atom_res != current_Resindex:
            current_Resindex = atom_res
            has_H, has_N = False, False


        if atom_res == current_Resindex:

            if atom.getName() == 'N':
                has_N = True
                N_coords = csx_obj.Vec_3D(atom.getCoords())

            elif atom.getName() == 'H':
                has_H = True
                H_coords = csx_obj.Vec_3D(atom.getCoords())

            if has_H and has_N:
                has_H, has_N = False, False
                vectors[atom_res] = csx_obj.Vec_3D(N_coords - H_coords).normalize()

    model_data.append(vectors)

S2_calced = {}

# az STR-ből származó S2 értékeken megy
for resnum in [int(s2rec.resnum) for s2rec in S2_records]:

    x2, y2, z2, xy, xz, yz = 0, 0, 0, 0, 0, 0

    # a modelleken megy
    for model in model_data:

        # adott modellben adott resnum-ra a normalizált vektorok koordinátái
        x, y, z = model[resnum].v[0], model[resnum].v[1], model[resnum].v[2]

        x2 += x ** 2
        y2 += y ** 2
        z2 += z ** 2
        xy += x * y
        xz += x * z
        yz += y * z

    x2 /= len(model_data)   # STR-ből az S2 adatok számával osztok
    y2 /= len(model_data)
    z2 /= len(model_data)
    xy /= len(model_data)
    xz /= len(model_data)
    yz /= len(model_data)

    s2 = 3 / 2.0 * (x2 ** 2 +
                    y2 ** 2 +
                    z2 ** 2 +
                    2 * xy ** 2 +
                    2 * xz ** 2 +
                    2 * yz ** 2) - 0.5

    S2_calced[resnum] = s2


# S2 correlation calc
def calcCorrel(S2_calced, S2_records):
    if len(S2_calced) != len(S2_records):
        return -2

    M = [0.0, 0.0, 0.0]
    D = [0.0, 0.0]

    for i, j in enumerate([int(s2rec.resnum) for s2rec in S2_records]):
        print(i)
        calc = S2_calced[j]
        exp  = float(S2_records[i].S2_value)

        M[0] += calc
        M[1] += exp
        M[2] += calc * exp

    M[0] /= len(S2_calced)
    M[1] /= len(S2_calced)
    M[2] /= len(S2_calced)

    for i, j in enumerate([int(s2rec.resnum) for s2rec in S2_records]):
        calc = S2_calced[j]
        exp  = float(S2_records[i].S2_value)

        D[0] += (calc - M[0]) ** 2
        D[1] += (exp  - M[1]) ** 2

    D[0] /= len(S2_calced)
    D[0] = math.sqrt(D[0])
    D[1] /= len(S2_calced)
    D[1] = math.sqrt(D[1])

    if D[0] * D[1] == 0:
        return -2
    else:
        return (M[2] - (M[0] * M[1])) / (D[0] * D[1])

print(calcCorrel(S2_calced, S2_records))
