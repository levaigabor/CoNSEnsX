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
import matplotlib.pyplot as plt
import nmrpystar

# own modules
import csx_libs.misc as misc
import csx_libs.methods as csx_func
from csx_libs.objects import *


version = "1.0"
pales   = "/home/daniel/Dokumente/önlab/gz_pack/pales/linux/pales"

parser = misc.createParser()
args = parser.parse_args()


if args.c:
    print(misc.cred)
    raise SystemExit

if args.help:
    print(misc.help_text)
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

#################
################# MODUL merge :D
#################

pales   = "/home/daniel/Dokumente/önlab/gz_pack/pales/linux/pales"
rdc_lc_model = "bic"

shortcodes = {
    'ALA':'A',  'ASP':'D',  'ASN':'N',  'ARG':'R',  'CYS':'C',  'GLY':'G',
    'GLU':'E',  'GLN':'Q',  'HIS':'H',  'ILE':'I',  'LEU':'L',  'LYS':'K',
    'MET':'M',  'PHE':'F',  'PRO':'P',  'SER':'S',  'THR':'T',  'TRP':'W',
    'TYR':'Y',  'VAL':'V'
}

def callPalesOn(pdb_files, RDC_list):
    try:
        os.remove("pales.out")                  # remove output file if present
    except OSError:
        pass

    for o, pdb_file in enumerate(pdb_files):
        #-------------------  Open file and read PDB data  -------------------#
        try:
            pdb_file  = "temp/" + pdb_file
            input_pdb = open(pdb_file)
        except IOError:
            print("Input file " + pdb_file + " was not found")
            raise SystemExit                    # exit if input file not found

        seg = []                                # list storing sequence data

        for line in input_pdb:
            if line.startswith("ATOM") and line.split()[2] == "CA":
                resname = line.split()[3]       # get residue name
                seg.append(resname)             # append new segname to list


        #-----------------------  Write sequence data  -----------------------#
        short_seg = ""

        for i in range(len(seg)):
            short_seg += shortcodes[seg[i]]

        my_line      = "DATA SEQUENCE "
        char_counter = 0
        row_counter  = 0
        pales_dummy  = open('pales_dummy.txt', 'w')

        for char in short_seg:
            if char_counter == 10:          # write aa output in 10 wide blocks
                my_line      += " "
                char_counter = 0
                row_counter  += 1

                if row_counter == 5:        # write 5 block per line
                    pales_dummy.write(my_line + "\n")
                    char_counter = 0
                    row_counter  = 0
                    my_line = "DATA SEQUENCE "

            my_line      += char
            char_counter += 1

        pales_dummy.write(my_line + "\n")    # write last line of aa output


        #-----------------------  Write dummy dipoles  -----------------------#
        pales_dummy.write(
        "\nVARS RESID_I RESNAME_I ATOMNAME_I " +
        "RESID_J RESNAME_J ATOMNAME_J D DD W\n" +
        "FORMAT %5d  %6s  %6s  %5d  %6s  %6s  %9.3f  %9.3f  %.2f \n\n"
        )

        for RDC_record in RDC_list:
            # print aligned dummy dipole output if present
            pales_dummy.write(
                "%5s  %6s  %6s  %5s  %6s  %6s  %9.3f  %9.3f  %.2f\n" % (
                str(RDC_record.resnum1) + 'A', seg[RDC_record.resnum1 - 1],
                str(RDC_record.atom1),
                str(RDC_record.resnum1) + 'A', seg[RDC_record.resnum1 - 1],
                str(RDC_record.atom2),
                RDC_record.RDC_value, 1.000,  1.00))

        pales_dummy.close()

        outfile = open("pales.out", 'a')    # open output file with append mode
        DEVNULL = open(os.devnull, 'w')     # open systems /dev/null

        #print("calculating: " + pdb_file)
        print("call Pales on model: " +
              str(o + 1) + '/' + str(len(pdb_files)), end="\r")
        sys.stdout.flush()


        if args.R:                          # if SVD is enabled
            subprocess.call([pales,
                            "-inD", "pales_dummy.txt",  # pales dummy file
                            "-pdb", pdb_file,           # pdb file
                            '-' + rdc_lc_model,         # rdc lc model
                            "-bestFit"],                # SVD
                            stdout = outfile,
                            stderr = DEVNULL)
        else:                               # if SVD is disabled (default)
            subprocess.call([pales,
                            "-inD", "pales_dummy.txt",  # pales dummy file
                            "-pdb", pdb_file,           # pdb file
                            '-' + rdc_lc_model],        # rdc lc model
                            stdout = outfile,
                            stderr = DEVNULL)
    print()


#####for RDC_list in RDC_lists:

callPalesOn(pdb_models, RDC_lists[0])

averageRDC = csx_func.avgPalesRDCs("pales.out")

for RDC_type in averageRDC.keys():
    print(RDC_type)




def calcCorrel(averageRDC, RDCtype, RDC_list):
    if len(averageRDC[RDCtype]) != len(RDC_lists[0]):
        return -2

    M = [0.0, 0.0, 0.0]
    D = [0.0, 0.0]

    for i, j in enumerate(averageRDC[RDCtype].keys()):
        calc = averageRDC[RDCtype][j]
        exp  = RDC_lists[0][i].RDC_value

        M[0] += calc
        M[1] += exp
        M[2] += calc * exp

    M[0] /= len(RDC_lists[0])
    M[1] /= len(RDC_lists[0])
    M[2] /= len(RDC_lists[0])

    for i, j in enumerate(averageRDC[RDCtype].keys()):
        calc = averageRDC[RDCtype][j]
        exp  = RDC_lists[0][i].RDC_value

        D[0] += (calc - M[0]) ** 2
        D[1] += (exp  - M[1]) ** 2

    D[0] /= len(RDC_lists[0])
    D[0] = math.sqrt(D[0])
    D[1] /= len(RDC_lists[0])
    D[1] = math.sqrt(D[1])

    if D[0] * D[1] == 0:
        return -2
    else:
        return (M[2] - (M[0] * M[1])) / (D[0] * D[1])


def calcQValue(averageRDC, RDCtype, RDC_list):

    if len(averageRDC[RDCtype]) != len(RDC_lists[0]):
        return -2

    D2, E2, C2 = 0, 0, 0

    for i, j in enumerate(averageRDC[RDCtype].keys()):
        calc = averageRDC[RDCtype][j]
        exp  = RDC_lists[0][i].RDC_value

        D2 += (calc - exp) ** 2
        E2 += exp ** 2

    Q = 100 * math.sqrt(D2) / math.sqrt(E2)

    return round(Q, 6)


def calcRMSD(averageRDC, RDCtype, RDC_list):

    if len(averageRDC[RDCtype]) != len(RDC_lists[0]):
        return -2

    D2 = 0

    for i, j in enumerate(averageRDC[RDCtype].keys()):
        calc = averageRDC[RDCtype][j]
        exp  = RDC_lists[0][i].RDC_value

        D2 += (calc - exp) ** 2

    RMSD = math.sqrt(D2 / len(RDC_lists[0]))

    return round(RMSD, 6)


def makeGraph(averageRDC, RDCtype, RDC_list):

    minrn = 0
    maxrn = len(averageRDC[RDCtype])

    min_calc = min(averageRDC[RDCtype].values())
    max_calc = max(averageRDC[RDCtype].values())

    exp_values = []
    for RDC_record in RDC_lists[0]:
        exp_values.append(RDC_record.RDC_value)

    min_exp = min(exp_values)
    max_exp = max(exp_values)

    miny = min(min_calc, min_exp)
    maxy = max(max_calc, max_exp)

    exp_line  = []
    calc_line = []

    print(min(averageRDC[RDCtype].keys()))
    print(max(averageRDC[RDCtype].keys()))


    for i, j in enumerate(averageRDC[RDCtype].keys()):
        calc = averageRDC[RDCtype][j]
        exp  = RDC_lists[0][i].RDC_value

        exp_line.append(exp)
        calc_line.append(calc)

    plt.plot(exp_line, linewidth=2.0, color='red', label='exp')
    plt.plot(calc_line, linewidth=2.0, color='blue', label='calc')
    plt.axis([minrn, maxrn, miny, maxy])
    plt.legend(loc='lower left')
    plt.show()


print("Correl: ", calcCorrel(averageRDC, '0_N_H', RDC_lists[0]))
print("Q-val:  ", calcQValue(averageRDC, '0_N_H', RDC_lists[0]))
print("RMSD:   ", calcRMSD(averageRDC, '0_N_H', RDC_lists[0]))
makeGraph(averageRDC, '0_N_H', RDC_lists[0])
