#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function
import os
import sys
import re
import subprocess
import math
import copy
import prody

import matplotlib.pyplot as plt
import numpy as np

from .objects import *


pales = "/home/daniel/Dokumente/önlab/gz_pack/pales/linux/pales"

shortcodes = {
    'ALA':'A',  'ASP':'D',  'ASN':'N',  'ARG':'R',  'CYS':'C',  'GLY':'G',
    'GLU':'E',  'GLN':'Q',  'HIS':'H',  'ILE':'I',  'LEU':'L',  'LYS':'K',
    'MET':'M',  'PHE':'F',  'PRO':'P',  'SER':'S',  'THR':'T',  'TRP':'W',
    'TYR':'Y',  'VAL':'V'
}


def pdb_cleaner(PDB_file):
    """
    Performs some basic formatting on the given PDB file to make it suitable
    for further calculations
    """
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
    """
    Split the given PDB file into models, each model becomes a separate
    PDB file placed in the "temp" folder
    """
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
    """
    Writes pales dummy from the given RDC values, and call Pales with the
    given parameters
    """
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
                RDC_record.value, 1.000,  1.00))

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
    """
    Returns a dictonary with the average RDCs:
    averageRDC[residue] = value
    """
    pales_out       = open(pales_out)
    n_of_structures = 0
    averageRDC      = {}
    npair           = 0

    for line in pales_out:
        if re.match("REMARK \d+ couplings", line):
            n_of_structures += 1                # n_of_structures to divide by

        elif re.match("\s+ \d+", line):
            resnum  = int(line.split()[0])
            resnum2 = int(line.split()[3])
            atom    = line.split()[2]
            atom2   = line.split()[5]
            D       = float(line.split()[8])    # D coloumn of pales output
            RDCtype = str(resnum2 - resnum) + "_" + atom + "_" + atom2

            if resnum in averageRDC.keys():
                averageRDC[resnum] += D
            else:
                averageRDC[resnum] = D

    for res_num in averageRDC.keys():
        averageRDC[res_num] /= n_of_structures

    return averageRDC


def calcS2(PDB_file, S2_records):
    """
    Returns a dictonary with the average S2 values:
    S2_calced[residue] = value
    """
    model_list = []
    model_num = 1

    while True:
        try:
            with suppress_output():
                # parsing PDB file into models (model_list)
                model_list.append(prody.parsePDB(PDB_file,
                                                 model=model_num, ter=True))
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
                    N_coords = Vec_3D(atom.getCoords())

                elif atom.getName() == 'H':
                    has_H = True
                    H_coords = Vec_3D(atom.getCoords())

                if has_H and has_N:
                    has_H, has_N = False, False
                    vectors[atom_res] = Vec_3D(N_coords -
                                               H_coords).normalize()

        model_data.append(vectors)

    S2_calced = {}

    # iterating over STR records
    for resnum in [int(s2rec.resnum) for s2rec in S2_records]:

        x2, y2, z2, xy, xz, yz = 0, 0, 0, 0, 0, 0

        # iterating over PDB models
        for m in model_data:

            # coordinates in model at a given resnum
            x, y, z = m[resnum].v[0], m[resnum].v[1], m[resnum].v[2]

            x2 += x ** 2
            y2 += y ** 2
            z2 += z ** 2
            xy += x * y
            xz += x * z
            yz += y * z

        x2 /= len(model_data)
        y2 /= len(model_data)
        z2 /= len(model_data)
        xy /= len(model_data)
        xz /= len(model_data)
        yz /= len(model_data)

        # S2 calcuation
        s2 = 3 / 2.0 * (x2 ** 2     + y2 ** 2     + z2 ** 2 +
                        2 * xy ** 2 + 2 * xz ** 2 + 2 * yz ** 2) - 0.5

        S2_calced[resnum] = s2

    return S2_calced


def calcCorrel(calced, experimental):
    """
    Calculates correlation between calculated and experimental data
    "calced" is a dict containing values for residues (as keys)
    "experimental" is a list containing STR record objects
    """
    if len(calced) != len(experimental):
        return -2

    M = [0.0, 0.0, 0.0]
    D = [0.0, 0.0]

    for i, j in enumerate(calced.keys()):
        calc = calced[j]
        exp  = experimental[i].value

        M[0] += calc
        M[1] += exp
        M[2] += calc * exp

    M[0] /= len(experimental)
    M[1] /= len(experimental)
    M[2] /= len(experimental)

    for i, j in enumerate(calced.keys()):
        calc = calced[j]
        exp  = experimental[i].value

        D[0] += (calc - M[0]) ** 2
        D[1] += (exp  - M[1]) ** 2

    D[0] /= len(experimental)
    D[0] = math.sqrt(D[0])
    D[1] /= len(experimental)
    D[1] = math.sqrt(D[1])

    if D[0] * D[1] == 0:
        return -2
    else:
        return (M[2] - (M[0] * M[1])) / (D[0] * D[1])


def calcQValue(calced, experimental):
    """
    Calculates Q-value for calculated and experimental data
    "calced" is a dict containing values for residues (as keys)
    "experimental" is a list containing STR record objects
    """
    if len(calced) != len(experimental):
        return -2

    D2, E2, C2 = 0, 0, 0

    for i, j in enumerate(calced.keys()):
        calc = calced[j]
        exp  = experimental[i].value

        D2 += (calc - exp) ** 2
        E2 += exp ** 2

    Q = 100 * math.sqrt(D2) / math.sqrt(E2)

    return round(Q, 6)


def calcRMSD(calced, experimental):
    """
    Calculates RMSD for calculated and experimental data
    "calced" is a dict containing values for residues (as keys)
    "experimental" is a list containing STR record objects
    """
    if len(calced) != len(experimental):
        return -2

    D2 = 0

    for i, j in enumerate(calced.keys()):
        calc = calced[j]
        exp  = experimental[i].value

        D2 += (calc - exp) ** 2

    RMSD = math.sqrt(D2 / len(experimental))

    return round(RMSD, 6)


def makeGraph(calced, my_experimental):
    """
    X axis -> residue numbers, Y axis -> values
    "calced" is a dict containing values for residues (as keys)
    "experimental" is a list containing STR record objects
    """
    experimental = copy.deepcopy(my_experimental)

    min_calc = min(calced.values())
    max_calc = max(calced.values())

    exp_values = []

    for record in experimental:
        exp_values.append(record.value)

    min_exp = min(exp_values)
    max_exp = max(exp_values)

    miny = min(min_calc, min_exp)               # get minimum value
    maxy = max(max_calc, max_exp)               # get maximum value

    exp_line, calc_line = [], []

    for k in range(0, max(calced.keys()) + 0):  # fetch data from arguments
        if k in calced.keys():
            calc = calced[k]
            exp  = experimental.pop(0).value

            exp_line.append(exp)
            calc_line.append(calc)

        else:
            exp_line.append(None)   # append 'None' where data is missing
            calc_line.append(None)

    # connect line over missing (None) values, more info at ->
    # http://stackoverflow.com/questions/14399689/
    # matplotlib-drawing-lines-between-points-ignoring-missing-data
    exp_line  = np.array(exp_line).astype(np.double)
    exp_mask  = np.isfinite(exp_line)
    calc_line = np.array(calc_line).astype(np.double)
    calc_mask = np.isfinite(calc_line)

    # x axis values as numpy array
    xs = np.arange(max(calced.keys())+2)
    # experimental values with 'None' values masked
    plt.plot(xs[exp_mask], exp_line[exp_mask],
             linewidth=2.0, color='red', marker='o', label='exp')
    # calculated values with 'None' values masked
    plt.plot(xs[calc_mask], calc_line[calc_mask],
             linewidth=2.0, color='blue', marker='o', label='calc')
    # setting axis limits
    plt.axis([min(calced.keys()) - 1, max(calced.keys()) + 1, miny, maxy])
    plt.legend(loc='lower left')
    plt.xlabel('residue number')
    plt.ylabel('value')
    #plt.show()


def makeCorrelGraph(calced, experimental):
    """
    X axis -> experimental values, Y axis -> calculated values
    "calced" is a dict containing values for residues (as keys)
    "experimental" is a list containing STR record objects
    """
    if len(calced) != len(experimental):
        return -2

    min_calc = min(calced.values())
    max_calc = max(calced.values())

    exp_values = []
    for record in experimental:
        exp_values.append(record.value)

    min_exp = min(exp_values)
    max_exp = max(exp_values)

    miny = min(min_calc, min_exp)               # get minimum value
    maxy = max(max_calc, max_exp)               # get maximum value

    exp_line, calc_line = [], []

    for i, j in enumerate(calced.keys()):       # fetch data from arguments
        calc = calced[j]
        exp  = experimental[i].value

        exp_line.append(exp)
        calc_line.append(calc)

    diag = []

    for i in np.arange(miny, maxy * 1.42, 0.1): # draw graph diagonal
        diag.append(i)

    plt.plot(diag, diag, linewidth=2.0, color='red')
    plt.plot(exp_line, calc_line, 'bo')
    plt.axis([miny, maxy, miny, maxy])
    plt.xlabel('experimental')
    plt.ylabel('calculated')
    #plt.show()
