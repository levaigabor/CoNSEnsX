#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import time
import os
import subprocess
import re
import math

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

parsed = nmrpystar.parse(myString)
if parsed.status == 'success':
    pass
    #print 'it worked!!  ', parsed.value
else:
    print 'Error during STR parsing: ', parsed


class RDC_Record(object):
    """Class for storing RDC data"""

    def __init__(self, resnum1, atom1, resnum2, atom2, RDC_value):
        self.RDC_type  = (str(int(resnum1) - int(resnum2))
                          + '_' + atom1 + '_' + atom2)
        self.resnum1   = int(resnum1)
        self.atom1     = atom1
        self.resnum2   = int(resnum2)
        self.atom2     = atom2
        self.RDC_value = float(RDC_value)


def get_RDC_lists(dataBlock):
    """Returns RDC lists as lists containing RDC_Record objects"""

    list_number = 1
    RDC_lists   = []

    while True:
        saveShiftName = 'RDC_list_' + str(list_number)
        try:
            saveShifts    = dataBlock.saves[saveShiftName]
        except KeyError:
            break
        loopShifts     = saveShifts.loops[-1]
        RDC_records    = []

        # STR key values recognised by this program
        rdc_types_keys = ["RDC.RDC_code", "Residual_dipolar_coupling_ID"]
        rdc_res1_keys  = ["RDC.Seq_ID_1", "Atom_one_residue_seq_code"]
        rdc_atom1_keys = ["RDC.Atom_type_1", "Atom_one_atom_name"]
        rdc_res2_keys  = ["RDC.Seq_ID_2", "Atom_two_residue_seq_code"]
        rdc_atom2_keys = ["RDC.Atom_type_2", "Atom_two_atom_name"]
        rdc_value_keys = ["RDC.Val", "Residual_dipolar_coupling_value"]

        for ix in range(len(loopShifts.rows)):   # fetch values from STR file
            row = loopShifts.getRowAsDict(ix)

            # for my_RDC_type in rdc_types_keys:   # fetch RDC type
            #     if my_RDC_type in row.keys():
            #         RDC_type = row[my_RDC_type]

            for my_resnum1 in rdc_res1_keys:     # fetch 1. residue number
                if my_resnum1 in row.keys():
                    resnum1 = row[my_resnum1]

            for my_atom1 in rdc_atom1_keys:      # fetch 1. atom name
                if my_atom1 in row.keys():
                    atom1 = row[my_atom1]

            for my_resnum2 in rdc_res2_keys:     # fetch 2. residue number
                if my_resnum2 in row.keys():
                    resnum2 = row[my_resnum2]

            for my_atom2 in rdc_atom2_keys:      # fetch 2. atom name
                if my_atom2 in row.keys():
                    atom2 = row[my_atom2]

            for my_RDC_value in rdc_value_keys:  # fetch RDC value
                if my_RDC_value in row.keys():
                    RDC_value = row[my_RDC_value]

            # check if all parameters are fetched
            if (resnum1 and atom1 and resnum2 and atom2 and RDC_value):
                # append RDC_Record object to list
                RDC_records.append(RDC_Record(resnum1, atom1,
                                              resnum2, atom2, RDC_value))
            else:
                raise ValueError('Unhandled value in STR file, check the\
                                  specifications')

        RDC_lists.append(RDC_records)
        list_number += 1

    return RDC_lists

RDC_lists  = get_RDC_lists(parsed.value)
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

    for pdb_file in pdb_files:
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


        # for i in range(len(seg)):
        #     if seg[i] == "PRO":         # skip PRO named aa-s
        #         continue

        #     # print aligned dummy dipole output
        #     pales_dummy.write(
        #         "%5s  %6s  %6s  %5s  %6s  %6s  %9.3f  %9.3f  %.2f\n" % (
        #         str(i+1)+'A', seg[i], 'H',
        #         str(i+1)+'A', seg[i], 'N',
        #         0.000, 1.000,  1.00))
        #     pales_dummy.write(
        #         "%5s  %6s  %6s  %5s  %6s  %6s  %9.3f  %9.3f  %.2f\n" % (
        #         str(i+1)+'A', seg[i], 'N',
        #         str(i+1)+'A', seg[i], 'C',
        #         0.000, 1.000,  1.00))
        #     pales_dummy.write(
        #         "%5s  %6s  %6s  %5s  %6s  %6s  %9.3f  %9.3f  %.2f\n" % (
        #         str(i+1)+'A', seg[i], 'C',
        #         str(i+1)+'A', seg[i], 'CA',
        #         0.000, 1.000,  1.00))

        pales_dummy.close()

        outfile = open("pales.out", 'a')    # open output file with append mode
        DEVNULL = open(os.devnull, 'w')     # open systems /dev/null

        print("calculating: " + pdb_file)

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


#####for RDC_list in RDC_lists:

#callPalesOn(pdb_models, RDC_lists[0])


def avgPalesRDCs(pales_out):
    pales_out       = open(pales_out)
    n_of_structures = 0
    averageRDC      = {}
    npair           = 0
    #resnum, exp, calc = [] # MIK EZEK??


    for line in pales_out:
        if re.match("REMARK \d+ couplings", line):
            n_of_structures += 1 # n_of_structures to divide by

        elif re.match("\s+ \d+", line):
            resnum1 = int(line.split()[0])
            resnum2 = int(line.split()[3])
            atom1   = line.split()[2]
            atom2   = line.split()[5]
            D       = float(line.split()[8])  # D coloumn of pales output
            RDCtype = str(resnum2 - resnum1) + "_" + atom1 + "_" + atom2

            if RDCtype in averageRDC.keys():
                if resnum1 in averageRDC[RDCtype].keys():
                    averageRDC[RDCtype][resnum1] += D
                else:
                    averageRDC[RDCtype][resnum1] = D
            else:
                averageRDC[RDCtype] = {}
                averageRDC[RDCtype][resnum1] = D



    for RDCtype in averageRDC.keys():
        for res_num in averageRDC[RDCtype].keys():
            averageRDC[RDCtype][res_num] /= n_of_structures


    #print averageRDC
    return averageRDC

averageRDC = avgPalesRDCs("pales.out")

# for RDC_type in averageRDC.keys():
#     print RDC_type


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


print calcCorrel(averageRDC, '0_N_H', RDC_lists[0])
print calcQValue(averageRDC, '0_N_H', RDC_lists[0])
print calcRMSD(averageRDC, '0_N_H', RDC_lists[0])
