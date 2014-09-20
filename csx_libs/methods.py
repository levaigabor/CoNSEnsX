#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function
import os
import sys
import re
import subprocess
import math

import matplotlib.pyplot as plt

from .objects import *

pales   = "/home/daniel/Dokumente/önlab/gz_pack/pales/linux/pales"

shortcodes = {
    'ALA':'A',  'ASP':'D',  'ASN':'N',  'ARG':'R',  'CYS':'C',  'GLY':'G',
    'GLU':'E',  'GLN':'Q',  'HIS':'H',  'ILE':'I',  'LEU':'L',  'LYS':'K',
    'MET':'M',  'PHE':'F',  'PRO':'P',  'SER':'S',  'THR':'T',  'TRP':'W',
    'TYR':'Y',  'VAL':'V'
}


def pdb_cleaner(PDB_file):
    try:
        input_pdb = open(PDB_file)
    except FileNotFoundError:
        print(PDB_file + " was not found")
        raise SystemExit

    my_pdb = open("my_pdb.pdb", 'w')

    for line in input_pdb:
        line = line.strip()
        line = re.sub('[+-] ', '  ', line)

        if line.startswith("ATOM"):

            name = line.split()[2].strip()

            if name is "Q": continue
            if name is "NH": name = "H"

            chars, numbers = [], []
            for i in name:
                try:
                    numbers.append(int(i))
                except ValueError:
                    chars.append(i)

            name = (''.join(str(i) for i in chars) +    # characters
                    ''.join(str(i) for i in numbers))   # numbers

            if   len(name) == 1: name = " " + name + "  "
            elif len(name) == 2: name = " " + name + " "
            elif len(name) == 3: name = " " + name

            my_pdb.write(line[:11] + " %4s" % name + line[16:21] +
                         'A' + line[22:] + "\n")
            continue

        my_pdb.write(line + "\n")

    input_pdb.close()
    my_pdb.close()

    os.remove(PDB_file)
    os.rename("my_pdb.pdb", PDB_file)


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
        rdc_res1_keys  = ["RDC.Seq_ID_1", "Atom_one_residue_seq_code"]
        rdc_atom1_keys = ["RDC.Atom_type_1", "Atom_one_atom_name"]
        rdc_res2_keys  = ["RDC.Seq_ID_2", "Atom_two_residue_seq_code"]
        rdc_atom2_keys = ["RDC.Atom_type_2", "Atom_two_atom_name"]
        rdc_value_keys = ["RDC.Val", "Residual_dipolar_coupling_value"]

        for ix in range(len(loopShifts.rows)):   # fetch values from STR file
            row = loopShifts.getRowAsDict(ix)

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


def callPalesOn(pdb_files, RDC_list, lc_model, SVD_enable):
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


        if SVD_enable:                          # if SVD is enabled
            subprocess.call([pales,
                            "-inD", "pales_dummy.txt",  # pales dummy file
                            "-pdb", pdb_file,           # pdb file
                            '-' + lc_model,             # rdc lc model
                            "-bestFit"],                # SVD
                            stdout = outfile,
                            stderr = DEVNULL)
        else:                               # if SVD is disabled (default)
            subprocess.call([pales,
                            "-inD", "pales_dummy.txt",  # pales dummy file
                            "-pdb", pdb_file,           # pdb file
                            '-' + lc_model],            # rdc lc model
                            stdout = outfile,
                            stderr = DEVNULL)
    print()


def avgPalesRDCs(pales_out):
    pales_out       = open(pales_out)
    n_of_structures = 0
    averageRDC      = {}
    npair           = 0

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

    return averageRDC


def calcCorrel(averageRDC, RDCtype, RDC_list):
    if len(averageRDC[RDCtype]) != len(RDC_list):
        return -2

    M = [0.0, 0.0, 0.0]
    D = [0.0, 0.0]

    for i, j in enumerate(averageRDC[RDCtype].keys()):
        calc = averageRDC[RDCtype][j]
        exp  = RDC_list[i].RDC_value

        M[0] += calc
        M[1] += exp
        M[2] += calc * exp

    M[0] /= len(RDC_list)
    M[1] /= len(RDC_list)
    M[2] /= len(RDC_list)

    for i, j in enumerate(averageRDC[RDCtype].keys()):
        calc = averageRDC[RDCtype][j]
        exp  = RDC_list[i].RDC_value

        D[0] += (calc - M[0]) ** 2
        D[1] += (exp  - M[1]) ** 2

    D[0] /= len(RDC_list)
    D[0] = math.sqrt(D[0])
    D[1] /= len(RDC_list)
    D[1] = math.sqrt(D[1])

    if D[0] * D[1] == 0:
        return -2
    else:
        return (M[2] - (M[0] * M[1])) / (D[0] * D[1])


def calcQValue(averageRDC, RDCtype, RDC_list):

    if len(averageRDC[RDCtype]) != len(RDC_list):
        return -2

    D2, E2, C2 = 0, 0, 0

    for i, j in enumerate(averageRDC[RDCtype].keys()):
        calc = averageRDC[RDCtype][j]
        exp  = RDC_list[i].RDC_value

        D2 += (calc - exp) ** 2
        E2 += exp ** 2

    Q = 100 * math.sqrt(D2) / math.sqrt(E2)

    return round(Q, 6)


def calcRMSD(averageRDC, RDCtype, RDC_list):

    if len(averageRDC[RDCtype]) != len(RDC_list):
        return -2

    D2 = 0

    for i, j in enumerate(averageRDC[RDCtype].keys()):
        calc = averageRDC[RDCtype][j]
        exp  = RDC_list[i].RDC_value

        D2 += (calc - exp) ** 2

    RMSD = math.sqrt(D2 / len(RDC_list))

    return round(RMSD, 6)


def makeGraph(averageRDC, RDCtype, RDC_list):

    minrn = 0
    maxrn = len(averageRDC[RDCtype])

    min_calc = min(averageRDC[RDCtype].values())
    max_calc = max(averageRDC[RDCtype].values())

    exp_values = []
    for RDC_record in RDC_list:
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
        exp  = RDC_list[i].RDC_value

        exp_line.append(exp)
        calc_line.append(calc)

    plt.plot(exp_line, linewidth=2.0, color='red', label='exp')
    plt.plot(calc_line, linewidth=2.0, color='blue', label='calc')
    plt.axis([minrn, maxrn, miny, maxy])
    plt.legend(loc='lower left')
    # plt.show()
