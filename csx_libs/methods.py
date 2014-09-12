#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import re

from .objects import *

#--------------------------  Bringing PDB to format  -------------------------#
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
