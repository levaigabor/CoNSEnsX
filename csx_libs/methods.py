import os
import sys
import re
import subprocess
import math
import copy
import string
import random
import prody
import time
import hashlib
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# installed modules
import nmrpystar
from .objects import *


shortcodes = {
    'ALA': 'A',  'ASP': 'D',  'ASN': 'N',  'ARG': 'R',  'CYS': 'C',  'GLY': 'G',
    'GLU': 'E',  'GLN': 'Q',  'HIS': 'H',  'ILE': 'I',  'LEU': 'L',  'LYS': 'K',
    'MET': 'M',  'PHE': 'F',  'PRO': 'P',  'SER': 'S',  'THR': 'T',  'TRP': 'W',
    'TYR': 'Y',  'VAL': 'V'
}

# Equation and coefficients from:
# Wang & Bax (1996) JACS 118:2483-2494. Table 1, NMR + X-ray data
Jcoup_dict1 = {
    'A': {"3JHNCB": 3.39,  "3JHNHA": 6.98,  "3JHNC": 4.32, "3JHAC": 3.75},
    'B': {"3JHNCB": -0.94, "3JHNHA": -1.38, "3JHNC": 0.84, "3JHAC": 2.19},
    'C': {"3JHNCB": 0.07,  "3JHNHA": 1.72,  "3JHNC": 0.00, "3JHAC": 1.28},
    'THETA': {
        "3JHNCB": math.radians(60), "3JHNHA": math.radians(-60),
        "3JHNC":  math.radians(0),  "3JHAC":  math.radians(-60)}  # RAD!
}

# J. Am. Chem. Soc., Vol. 119, No. 27, 1997; Table 2 -> solution
Jcoup_dict2 = {
    'A': {"3JHNCB": 3.06,  "3JHNHA": 7.13, "3JHNC": 4.19, "3JHAC": 3.84},
    'B': {"3JHNCB": -0.74, "3JHNHA": 1.31, "3JHNC": 0.99, "3JHAC": 2.19},
    'C': {"3JHNCB": 0.10,  "3JHNHA": 1.56, "3JHNC": 0.03, "3JHAC": 1.20},
    'THETA': {
        "3JHNCB": math.radians(60),  "3JHNHA": math.radians(-60),
        "3JHNC":  math.radians(180), "3JHAC":  math.radians(120)}  # RAD!
}
# https://x86.cs.duke.edu/~brd/Teaching/Bio/asmb/Papers/NMR/nilges-jmr05.pdf
Jcoup_dict3 = {
    'A': {"3JHNCB": 3.26,  "3JHNHA": 7.13,  "3JHNC": 4.19, "3JHAC": 3.84},
    'B': {"3JHNCB": -0.87, "3JHNHA": -1.31, "3JHNC": 0.99, "3JHAC": 2.19},
    'C': {"3JHNCB": 0.10,  "3JHNHA": 1.56,  "3JHNC": 0.03, "3JHAC": 1.20},
    'THETA': {
        "3JHNCB": math.radians(60), "3JHNHA": math.radians(-60),
        "3JHNC":  math.radians(0),  "3JHAC":  math.radians(-60)}  # RAD!
}


def timeit(method):
    """Timer decorator to keep an eye on CPU hungry processes"""
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        print('\x1b[31m%r -> %2.2f sec\x1b[0m' % (method.__name__, te-ts),
              file=sys.stderr)
        return result

    return timed


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def getID(args):
    """Generates unique ID for calculation"""
    chars = string.ascii_uppercase + string.digits

    if args.i:
        my_id = args.i
    else:
        if os.path.exists("calculations"):
            while True:
                my_id    = ''.join(random.choice(chars) for _ in range(6))
                used_IDs = os.listdir("calculations")
                if my_id not in used_IDs:
                    break
        else:
            my_id = ''.join(random.choice(chars) for _ in range(6))
    my_path = "calculations/" + my_id + '/'
    print("Job started with ID: \033[0;35m" + my_id + "\033[0m")
    return my_id, my_path


def check_3rd_party(install_dir):
    config_file_name = install_dir + "/.config"
    if os.path.isfile(config_file_name):
        ThirdParty.get_thirdparty(config_file_name)
    else:
        init_conf = (
            "# CoNSEnsX config file\n" +
            "# Please provide full paths\n" +
            "pales=''\n" + "shiftx=''\n" +
            "prideDB=''\n" + "prideNMR=''"
        )

        init_conf_file = open(".config", "w")
        init_conf_file.write(init_conf)
        print("Please edit '.config' in the CoNSEnsX install directory")
        raise SystemExit


def get_PDB(args):
    """Gets PDB file or downloads PDF file from rcsb.org"""
    if args.PDB_file:
        my_PDB = args.PDB_file
    else:
        my_PDB = prody.fetchPDB(args.PDB_fetch, compressed=False)
        print()
    return my_PDB


@timeit
def get_model_list(PDB_file, PDB_file_name=None, ):
    """Parsing PDB file into models in the PDB_model object"""
    model_num = 1

    my_PDB_hash = hashlib.md5(open(PDB_file_name, 'rb').read()).hexdigest()

    try:
        prev_PDB_hash = open(".prevPDBhash").read().strip()
        if my_PDB_hash == prev_PDB_hash:
            print("PDB file is the same as in the previous run")
            print("Skip parsing PDB...")
            pickled_models = pickle.load(open("PDB_model.pickle", 'rb'))
            PDB_model.model_list = pickled_models
            return
    except FileNotFoundError:
        pass

    with open(".prevPDBhash", "w") as prev_HASH:
        prev_HASH.write(my_PDB_hash)

    while True:
        try:
            PDB_model(prody.parsePDB(PDB_file, model=model_num, ter=True))
            model_num += 1
        except prody.proteins.pdbfile.PDBParseError:
            break
        finally:
            pickle.dump( PDB_model.model_list, open("PDB_model.pickle", "wb"))


@timeit
def pdb_cleaner(PDB_file):
    """Performs some basic formatting on the given PDB file to make it suitable
       for further calculations (from RCSB to BMRB PDF format)"""
    try:
        input_pdb = open(PDB_file)
    except FileNotFoundError:
        print(PDB_file + " was not found")
        raise SystemExit

    my_pdb = open("my_pdb.pdb", 'w')
    max_resnum = 0
    min_resnum = 100000

    for line in input_pdb:
        line = line.strip()
        # line = re.sub('[+-] ', '  ', line)

        if line.startswith("ATOM"):

            name   = line.split()[2].strip()
            try:
                resnum = int(line.split()[5].strip())
            except ValueError:
                resnum = int(line.split()[4].strip())

            if resnum >= max_resnum:
                max_resnum = resnum
            if resnum < min_resnum:
                min_resnum = resnum

            if name == "Q":
                continue
            if name == "NH":
                name = "H"
            if name == "HN":
                name = "H"

            chars, numbers = [], []
            for i in name:
                try:
                    numbers.append(int(i))
                except ValueError:
                    chars.append(i)

            name = (''.join(str(i) for i in chars) +    # characters
                    ''.join(str(i) for i in numbers))   # numbers

            if len(name) == 1:
                name = " " + name + "  "
            elif len(name) == 2:
                name = " " + name + " "
            elif len(name) == 3:
                name = " " + name

            my_pdb.write(line[:11] + " %4s" % name + line[16:21] +
                         'A' + line[22:] + "\n")
            continue

        my_pdb.write(line + "\n")

    input_pdb.close()
    my_pdb.close()

    os.remove(PDB_file)
    os.rename("my_pdb.pdb", PDB_file)

    # set max resnum value for CSV_buffer class
    CSV_buffer.max_resnum = max_resnum
    CSV_buffer.min_resnum = min_resnum


def pdb_splitter(my_path, PDB_file):
    """Split the given PDB file into models, each model becomes a separate
       PDB file placed in the "temp" folder"""
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
            # replace oxygen names in line
            try:
                if line.split()[2] == "OC1":
                    line = line.replace('OC1', 'O  ')
                elif line.split()[2] == "OC2":
                    line = line.replace('OC2', 'OXT')
            except IndexError:
                pass
            my_data.append(line.strip())
        elif line.startswith("ENDMDL"):
            model_names.append(my_name)
            model_data.append(my_data)
            my_name = ""
            my_data = []
        else:
            continue

    my_pdb.close()

    for i in range(len(model_names)):
        file_name = my_path + "/model_" + model_names[i] + ".pdb"
        temp_pdb  = open(file_name, 'w')
        temp_pdb.write("HEADER    MODEL " + model_names[i] + "\n")

        for _ in model_data[i]:
            temp_pdb.write(_ + "\n")

        temp_pdb.write("END")
        temp_pdb.close()

    return len(model_names)


@timeit
def getNOE(NOE_file):
    NOE_key   = "_Gen_dist_constraint_list.Constraint_type"
    NOE_value = "NOE"
    in_frame  = False
    in_loop   = False
    loops     = []
    NOE_data  = []

    for line in open(NOE_file):
        tok = line.strip().split()

        if not tok:
            continue

        if tok[0] == NOE_key and tok[1] == NOE_value:
            in_frame = True
            continue

        if in_frame and tok[0] == "save_":
            break

        if tok[0] == "loop_":
            loop_keys = []
            loop_data = []
            in_loop = True
            continue

        if in_loop is True and tok[0] == "stop_":
            in_loop = False
            loops.append([loop_keys, loop_data])

        if in_loop:
            if tok[0].startswith('_'):
                loop_keys.append(tok[0])
            else:
                loop_data.append(tok)

    for loop in loops:
        try:
            ind_ID    = loop[0].index("_Gen_dist_constraint.ID")
            ind_seg1  = loop[0].index("_Gen_dist_constraint.Seq_ID_1")
            ind_seg2  = loop[0].index("_Gen_dist_constraint.Seq_ID_2")
            ind_comp1 = loop[0].index("_Gen_dist_constraint.Comp_ID_1")
            ind_comp2 = loop[0].index("_Gen_dist_constraint.Comp_ID_2")
            ind_atom1 = loop[0].index("_Gen_dist_constraint.Atom_ID_1")
            ind_atom2 = loop[0].index("_Gen_dist_constraint.Atom_ID_2")
            ind_bnd   = loop[0].index(
                "_Gen_dist_constraint.Distance_upper_bound_val")
        except ValueError:
            continue

        for data in loop[1]:
            # print(data[ind_ID], data[ind_seg1], data[ind_seg2],
            #       data[ind_comp1], data[ind_comp2], data[ind_atom1],
            #       data[ind_atom2], data[ind_bnd])
            NOE_data.append([data[ind_ID],    data[ind_seg1],  data[ind_seg2],
                             data[ind_comp1], data[ind_comp2], data[ind_atom1],
                             data[ind_atom2], data[ind_bnd]])

    return NOE_data


@timeit
def parseSTR(STR_file):
    """Parse BMRB file into a python object"""

    my_STR_hash = hashlib.md5(open(STR_file, 'rb').read()).hexdigest()

    try:
        prev_STR_hash = open(".prevSTRhash").read().strip()

        if my_STR_hash == prev_STR_hash:
            print("STR file is the same as in the previous run")
            print("Skip parsing STR...")
            parsed = pickle.load(open(".prevSTR.pickle", 'rb'))
            return parsed

    except FileNotFoundError:
        pass

    with open(".prevSTRhash", "w") as prev_HASH:
        prev_HASH.write(my_STR_hash)

    star_file = open(STR_file)         # open STR file
    myString = ""

    for line in star_file:                  # rean STR file into a string
        myString += line

    star_file.close()
    parsed = nmrpystar.parse(myString)      # parsing -> parsed.value

    if parsed.status != 'success':          # check if parsing was successful
        print('Error during STR parsing: ', parsed)
        raise SystemExit
    else:
        pickle.dump( parsed, open(".prevSTR.pickle", "wb"))
        return parsed


def get_RDC_lists(parsed_value):
    """Returns RDC lists as dictonaries containing RDC_Record objects,
       grouped by RDCtype (keys())"""
    list_number = 1
    RDC_lists   = []

    while True:
        saveShiftName = 'RDC_list_' + str(list_number)
        try:
            saveShifts = parsed_value.saves[saveShiftName]
        except KeyError:
            break
        loopShifts  = saveShifts.loops[-1]
        RDC_records = []

        # STR key values recognised by this program
        rdc_res1_keys  = ["RDC.Seq_ID_1", "Atom_one_residue_seq_code"]
        rdc_atom1_keys = ["RDC.Atom_type_1", "Atom_one_atom_name"]
        rdc_res2_keys  = ["RDC.Seq_ID_2", "Atom_two_residue_seq_code"]
        rdc_atom2_keys = ["RDC.Atom_type_2", "Atom_two_atom_name"]
        rdc_value_keys = ["RDC.Val", "Residual_dipolar_coupling_value"]

        for ix in range(len(loopShifts.rows)):   # fetch values from STR file
            row = loopShifts.getRowAsDict(ix)

            for my_resnum1 in rdc_res1_keys:     # fetch 1. residue number
                if my_resnum1 in list(row.keys()):
                    resnum1 = row[my_resnum1]

            for my_atom1 in rdc_atom1_keys:      # fetch 1. atom name
                if my_atom1 in list(row.keys()):
                    atom1 = row[my_atom1]

            for my_resnum2 in rdc_res2_keys:     # fetch 2. residue number
                if my_resnum2 in list(row.keys()):
                    resnum2 = row[my_resnum2]

            for my_atom2 in rdc_atom2_keys:      # fetch 2. atom name
                if my_atom2 in list(row.keys()):
                    atom2 = row[my_atom2]

            for my_RDC_value in rdc_value_keys:  # fetch RDC value
                if my_RDC_value in list(row.keys()):
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

    # split list into dict according to RDC types
    new_RDC_list = []
    for list_num, RDC_list in enumerate(RDC_lists):
        prev_type = ""
        RDC_dict  = {}

        for record in RDC_list:
            if prev_type != record.RDC_type:
                RDC_dict[record.RDC_type] = []
                RDC_dict[record.RDC_type].append(record)
            else:
                RDC_dict[record.RDC_type].append(record)

            prev_type = record.RDC_type

        new_RDC_list.append(RDC_dict)

    return new_RDC_list


def parseS2_STR(parsed_value):
    """Returns a dictonary with the parsed S2 data"""
    try:
        saveShifts = parsed_value.saves["order_param"]

        loopShifts = saveShifts.loops[-1]
        S2_records = []

        for ix in range(len(loopShifts.rows)):   # fetch values from STR file
            row = loopShifts.getRowAsDict(ix)

            S2_records.append(S2_Record(row["Residue_seq_code"],
                                        row["Atom_name"],
                                        row["S2_value"]))

        # split list into dict according to S2 types
        S2_dict = {}
        prev_type = ""

        for record in S2_records:
            if prev_type != record.type:
                S2_dict[record.type] = []
                S2_dict[record.type].append(record)
            else:
                S2_dict[record.type].append(record)

            prev_type = record.type

        return S2_dict

    except KeyError:
        print("No S2 parameter list found")
        return None


def parse_sidechain_S2_STR(parsed_value):
    """Returns a dictonary with the parsed S2 data"""
    try:
        saveShifts = parsed_value.saves["side-chain_methyl_order_parameters"]

        loopShifts = saveShifts.loops[-1]
        S2_records = []

        for ix in range(len(loopShifts.rows)):   # fetch values from STR file
            row = loopShifts.getRowAsDict(ix)

            S2_records.append(S2_Record(row["Residue_seq_code"],
                                        row["Atom_name"],
                                        row["S2_value"]))

        return S2_records

    except KeyError:
        return None


def parseJcoup_STR(parsed_value):
    """Returns a dictonary with the parsed J-coupling data"""
    try:
        saveShifts = parsed_value.saves["coupling_constants"]

        loopShifts = saveShifts.loops[-1]
        jcoup_records = []

        for ix in range(len(loopShifts.rows)):   # fetch values from STR file
            row = loopShifts.getRowAsDict(ix)

            jcoup_records.append(JCoup_Record(row["Atom_one_residue_seq_code"],
                                              row["Coupling_constant_code"],
                                              row["Coupling_constant_value"]))

        # split list into dict according to J-cuopling types
        jcoup_dict = {}
        prev_type = ""

        for record in jcoup_records:
            if prev_type != record.type:
                jcoup_dict[record.type] = []
                jcoup_dict[record.type].append(record)
            else:
                jcoup_dict[record.type].append(record)

            prev_type = record.type

        return jcoup_dict

    except KeyError:
        print("No J-coupling parameter list found")
        return None


def parseChemShift_STR(parsed_value):
    """Returns ChemShift lists as dictonaries containing ChemShift_Record
    objects, grouped by Atom_name (keys())"""
    list_number = 1
    ChemShift_lists = []

    while True:
        saveShiftName = 'chem_shift_list_' + str(list_number)
        try:
            saveShifts = parsed_value.saves[saveShiftName]
        except KeyError:
            break

        loopShifts = saveShifts.loops[-1]
        ChemShift_records = []
        HA_sum = 0.0

        for ix in range(len(loopShifts.rows)):   # fetch values from STR file
            row = loopShifts.getRowAsDict(ix)

            if row["Atom_name"] == "HA2":
                HA_sum += float(row["Chem_shift_value"])
                continue

            if row["Atom_name"] == "HA3":
                HA_sum += float(row["Chem_shift_value"])

                ChemShift_records.append(
                    ChemShift_Record(
                        row["Residue_seq_code"], row["Residue_label"],
                        "HA", HA_sum / 2)
                )
                HA_sum = 0.0
                continue

            if row["Atom_name"] in ["HA", "CA", "N", "H"]:
                ChemShift_records.append(
                    ChemShift_Record(
                        row["Residue_seq_code"],
                        row["Residue_label"],
                        row["Atom_name"],
                        row["Chem_shift_value"])
                )

        ChemShift_lists.append(ChemShift_records)
        list_number += 1

    new_CS_list = []

    for ChemShift_list in ChemShift_lists:
        ChemShift_dict = {}

        for record in ChemShift_list:

            if record.atom_name in list(ChemShift_dict.keys()):
                ChemShift_dict[record.atom_name].append(record)
            else:
                ChemShift_dict[record.atom_name] = []
                ChemShift_dict[record.atom_name].append(record)

        new_CS_list.append(ChemShift_dict)

    return new_CS_list


@timeit
def callPalesOn(my_path, pdb_files, RDC_dict, lc_model, SVD_enable):
    """Writes pales dummy from the given RDC values, and call Pales with the
    given parameters"""
    for o, pdb_file in enumerate(pdb_files):
        #-------------------  Open file and read PDB data  -------------------#
        try:
            pdb_file  = my_path + pdb_file
            input_pdb = open(pdb_file)
        except IOError:
            print("Input file " + pdb_file + " was not found")
            raise SystemExit                    # exit if input file not found

        seg = []                                # list storing sequence data

        for line in input_pdb:
            if line.startswith("ATOM") and line.split()[2] == "CA":
                resname = line.split()[3]       # get residue name
                seg.append(resname)             # append new segname to list

        input_pdb.close()

        #-----------------------  Write sequence data  -----------------------#
        short_seg = ""

        for i in range(len(seg)):
            short_seg += shortcodes[seg[i]]

        my_line      = "DATA SEQUENCE "
        char_counter = 0
        row_counter  = 0
        pales_dummy  = open(my_path + 'pales_dummy.txt', 'w')

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

        lists = []
        for RDC_list in list(RDC_dict.keys()):
            lists.append(RDC_dict[RDC_list])

        for RDC_set in lists:
            for RDC_record in RDC_set:

                pales_dummy.write(
                    "%5s  %6s  %6s  %5s  %6s  %6s  %9.3f  %9.3f  %.2f\n" % (
                    str(RDC_record.resnum) + 'A', seg[RDC_record.resnum - 1],
                    str(RDC_record.atom),
                    str(RDC_record.resnum2) + 'A', seg[RDC_record.resnum2 - 1],
                    str(RDC_record.atom2),
                    RDC_record.value, 1.000,  1.00))

        pales_dummy.close()

        outfile = open(my_path + "pales.out", 'a')
        DEVNULL = open(os.devnull, 'w')

        print("call Pales on model: " +
              str(o + 1) + '/' + str(len(pdb_files)), end="\r")
        sys.stdout.flush()

        if SVD_enable:                          # if SVD is enabled
            subprocess.call([ThirdParty.pales,
                            "-inD", my_path + "pales_dummy.txt",
                            "-pdb", pdb_file,           # pdb file
                            '-' + lc_model,             # rdc lc model
                            "-bestFit"],                # SVD
                            stdout=outfile,
                            stderr=DEVNULL)
        else:                               # if SVD is disabled (default)
            subprocess.call([ThirdParty.pales,
                            "-inD", my_path + "pales_dummy.txt",
                            "-pdb", pdb_file,           # pdb file
                            '-' + lc_model],            # rdc lc model
                            stdout=outfile,
                            stderr=DEVNULL)
        outfile.close()
        DEVNULL.close()

    print()


@timeit
def callShiftxOn(my_path, pdb_files):
    """Call ShiftX on PDB models. Each output is appended to 'out_name'"""
    for i, pdb_file in enumerate(pdb_files):
        pdb_file = my_path + pdb_file
        out_name = my_path + "/modell_" + str(i+1) + ".out"
        subprocess.call([ThirdParty.shiftx, '1', pdb_file, out_name])

    averageHA, averageH, averageN, averageCA = {}, {}, {}, {}
    modHA, modH, modN, modCA                 = {}, {}, {}, {}
    model_data_list = []

    for some_file in os.listdir(my_path):
        if some_file.startswith("modell") and some_file.endswith(".out"):
            out_file = open(my_path + some_file)
            part = 0

            for line in out_file:
                if line.strip().startswith("NUM"):
                    part += 1

                    if part == 2:
                        break

                    continue

                if line.strip().startswith("---"):
                    continue

                if line.strip():
                    line_values = line.split()
                    resnum = int(line_values[0])
                    HA = float(line_values[2])
                    H  = float(line_values[3])
                    N  = float(line_values[4])
                    CA = float(line_values[5])

                    modHA[resnum] = HA
                    modH[resnum]  = H
                    modN[resnum]  = N
                    modCA[resnum] = CA

                    if resnum in list(averageHA.keys()):
                        averageHA[resnum] += HA
                    else:
                        averageHA[resnum] = HA

                    if resnum in list(averageH.keys()):
                        averageH[resnum] += H
                    else:
                        averageH[resnum] = H

                    if resnum in list(averageN.keys()):
                        averageN[resnum] += N
                    else:
                        averageN[resnum] = N

                    if resnum in list(averageCA.keys()):
                        averageCA[resnum] += CA
                    else:
                        averageCA[resnum] = CA

            out_file.close()
            model_data_list.append({"HA": modHA, "H": modH,
                                    "N":  modN, "CA": modCA})
            modHA, modH, modN, modCA = {}, {}, {}, {}

    for avg_dict in [averageHA, averageH, averageN, averageCA]:
        for key in avg_dict:
            avg_dict[key] /= len(pdb_files)

    return {"HA": averageHA, "H":  averageH,
            "N":  averageN,  "CA": averageCA}, model_data_list


def avgPalesRDCs(pales_out, my_RDC_type):
    """Returns a dictonary with the average RDCs for a given RDC type:
       averageRDC[residue] = value
       and calculated model data as a list of dictonaries
       model_data_list[{1: value}, ...]"""
    pales_out       = open(pales_out)
    n_of_structures = 0
    averageRDC      = {}
    model_data_list = []
    model_data_dict = {}
    first_run       = True

    for line in pales_out:
        if re.match("REMARK \d+ couplings", line):
            if first_run:
                first_run = False
                continue

            n_of_structures += 1                # n_of_structures to divide by

            model_data_list.append(model_data_dict)

            if model_data_dict:
                model_data_dict = {}

        elif re.match("\s+ \d+", line):
            resnum  = int(line.split()[0])
            resnum2 = int(line.split()[3])
            atom    = line.split()[2]
            atom2   = line.split()[5]
            D       = float(line.split()[8])    # D coloumn of pales output
            RDCtype = str(resnum2 - resnum) + "_" + atom + "_" + atom2

            # skip non relevant RDC data in the pales output file
            if my_RDC_type != RDCtype:
                continue

            if resnum in list(averageRDC.keys()):
                averageRDC[resnum] += D
            else:
                averageRDC[resnum] = D

            model_data_dict[resnum] = D

    model_data_list.append(model_data_dict)
    n_of_structures += 1
    pales_out.close()

    for res_num in list(averageRDC.keys()):
        averageRDC[res_num] /= n_of_structures

    return averageRDC, model_data_list


@timeit
def calcS2(model_list, S2_records, S2_type, fit, fit_range):
    """Returns a dictonary with the average S2 values:
    S2_calced[residue] = value"""

    # fitting models
    reference = model_list[0]

    if fit and not PDB_model.is_fitted:
        print("Start FITTING")
        for i in range(1, len(model_list)):
            mobile    = model_list[i]
            matches   = prody.matchChains(reference, mobile)
            match     = matches[0]
            ref_chain = match[0]
            mob_chain = match[1]

            if fit_range:
                weights = np.zeros((len(ref_chain), 1), dtype=np.int)

                fit_start, fit_end = fit_range.split('-')

                for i in range(int(fit_start) - 1, int(fit_end) - 1):
                    weights[i] = 1

            else:
                weights = np.ones((len(ref_chain), 1), dtype=np.int)

            t = prody.calcTransformation(mob_chain, ref_chain, weights)
            t.apply(mobile)

        PDB_model.is_fitted = True

    # get NH vectors from models (model_data[] -> vectors{resnum : vector})
    model_data = []
    s2_pairs   = {'N':  'H',
                  'CA': 'HA'}

    for model in model_list:
        current_Resindex = 1
        has_first, has_second = False, False
        vectors = {}

        for atom in model:
            # why not .getResnum() ???
            atom_res = atom.getResindex() + 1

            if atom_res != current_Resindex:
                current_Resindex = atom_res
                has_first, has_second = False, False

            if atom_res == current_Resindex:
                if atom.getName() == S2_type:
                    has_second = True
                    N_coords = Vec_3D(atom.getCoords())

                elif atom.getName() == s2_pairs[S2_type]:
                    has_first = True
                    H_coords = Vec_3D(atom.getCoords())

                if has_first and has_second:
                    has_first, has_second = False, False
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


@timeit
def calcDihedAngles(model_list):
    """Calculates backbone diherdral angles
       note: all returned angle values are in radian"""
    model_list = PDB_model.model_list

    JCoup_dicts = []

    for model_num, model in enumerate(model_list):
        current_Resindex = 1
        prev_C, my_N, my_CA, my_C = None, None, None, None
        JCoup_dict = {}

        for atom in model:
            atom_res = atom.getResindex() + 1

            if atom_res != current_Resindex:

                if (prev_C is not None and my_N is not None and
                        my_CA is not None and my_C is not None):

                    NCA_vec = my_N - my_CA
                    CN_vec  = prev_C - my_N
                    CCA_vec = my_C - my_CA

                    first_cross  = Vec_3D.cross(CN_vec, NCA_vec)
                    second_cross = Vec_3D.cross(CCA_vec, NCA_vec)

                    angle = Vec_3D.dihedAngle(first_cross, second_cross)

                    # reference for setting sign of angle
                    reference = Vec_3D.cross(first_cross, second_cross)

                    r1 = reference.normalize()
                    r2 = NCA_vec.normalize()

                    if ((r1 - r2).magnitude() < r2.magnitude()):
                        angle *= -1

                    JCoup_dict[current_Resindex] = -1 * math.radians(angle)

                current_Resindex = atom_res
                prev_C = my_C
                my_N, my_CA, my_C = None, None, None

            if atom_res == current_Resindex:
                if atom.getName() == 'N':
                    my_N = Vec_3D(atom.getCoords())
                elif atom.getName() == 'CA':
                    my_CA = Vec_3D(atom.getCoords())
                elif atom.getName() == 'C':
                    my_C = Vec_3D(atom.getCoords())

        JCoup_dicts.append(JCoup_dict)

    return JCoup_dicts


def calcPeptideBonds(PDB_file):
    """Calculates backbone diherdral angles (OMEGA) CA-N-C'-CA"""
    model_list = PDB_model.model_list
    dihedral_angles = {"<2":    0, "2-5": 0, "5-10": 0,
                       "10-20": 0, ">20": 0}

    for model_num, model in enumerate(model_list):
        current_Resindex = 1
        prev_CA, my_N, my_CA, my_C = None, None, None, None

        for atom in model:
            atom_res = atom.getResindex() + 1

            if atom_res != current_Resindex:

                if (prev_CA is not None and my_N is not None and
                        my_CA is not None and prev_C is not None):

                    NCA_vec = my_N - my_CA
                    CN_vec  = prev_CA - my_N
                    CCA_vec = prev_C - my_CA

                    first_cross  = Vec_3D.cross(CN_vec, NCA_vec)
                    second_cross = Vec_3D.cross(CCA_vec, NCA_vec)

                    angle = Vec_3D.dihedAngle(first_cross, second_cross)

                    # reference for setting sign of angle
                    reference = Vec_3D.cross(first_cross, second_cross)

                    r1 = reference.normalize()
                    r2 = NCA_vec.normalize()

                    if ((r1 - r2).magnitude() < r2.magnitude()):
                        angle *= -1

                    if abs(angle) < 2:
                        dihedral_angles["<2"] += 1
                    elif abs(angle) < 5:
                        dihedral_angles["2-5"] += 1
                    elif abs(angle) < 10:
                        dihedral_angles["5-10"] += 1
                    elif abs(angle) < 20:
                        dihedral_angles["10-20"] += 1
                    else:
                        dihedral_angles[">20"] += 1

                current_Resindex = atom_res
                prev_CA = my_CA
                prev_C =  my_C
                my_N, my_CA, my_C = None, None, None

            if atom_res == current_Resindex:
                if atom.getName() == 'N':
                    my_N = Vec_3D(atom.getCoords())
                elif atom.getName() == 'CA':
                    my_CA = Vec_3D(atom.getCoords())
                elif atom.getName() == 'C':
                    my_C = Vec_3D(atom.getCoords())

    print("Peptide (CA-N-C'-CA) bond angle distribution:")
    print("   <2 -> " + str(dihedral_angles["<2"]))
    print("  2-5 -> " + str(dihedral_angles["2-5"]))
    print(" 5-10 -> " + str(dihedral_angles["5-10"]))
    print("10-20 -> " + str(dihedral_angles["10-20"]))
    print("  >20 -> " + str(dihedral_angles[">20"]))


def calcNH_Angles(PDB_file):
    """Calculates backbone diherdral angles (OMEGA) H-N-C=O"""
    model_list = PDB_model.model_list
    dihedral_angles = {"<2":    0, "2-5": 0, "5-10": 0,
                       "10-20": 0, ">20": 0}

    for model_num, model in enumerate(model_list):
        current_Resindex = 1
        prev_O, prev_C, my_N, my_H = None, None, None, None

        for atom in model:
            atom_res = atom.getResindex() + 1

            if atom_res != current_Resindex:

                if (prev_O is not None and prev_C is not None and
                        my_N is not None and my_H is not None):

                    NH_vec = my_H - my_N
                    CN_vec = my_N - prev_C
                    OC_vec = prev_C - prev_O

                    first_cross  = Vec_3D.cross(NH_vec, CN_vec)
                    second_cross = Vec_3D.cross(OC_vec, CN_vec)

                    angle = Vec_3D.dihedAngle(first_cross, second_cross)

                    # reference for setting sign of angle
                    reference = Vec_3D.cross(first_cross, second_cross)

                    r1 = reference.normalize()
                    r2 = NH_vec.normalize()

                    if ((r1 - r2).magnitude() < r2.magnitude()):
                        angle *= -1

                    if abs(angle) < 2:
                        dihedral_angles["<2"] += 1
                    elif abs(angle) < 5:
                        dihedral_angles["2-5"] += 1
                    elif abs(angle) < 10:
                        dihedral_angles["5-10"] += 1
                    elif abs(angle) < 20:
                        dihedral_angles["10-20"] += 1
                    else:
                        dihedral_angles[">20"] += 1

                current_Resindex = atom_res
                prev_O = my_O
                prev_C = my_C
                my_N, my_H = None, None

            if atom_res == current_Resindex:
                if atom.getName() == 'N':
                    my_N = Vec_3D(atom.getCoords())
                elif atom.getName() == 'H':
                    my_H = Vec_3D(atom.getCoords())
                elif atom.getName() == 'C':
                    my_C = Vec_3D(atom.getCoords())
                elif atom.getName() == 'O':
                    my_O = Vec_3D(atom.getCoords())

    print("Peptide (H-N-C=O) bond angle distribution:")
    print("   <2 -> " + str(dihedral_angles["<2"]))
    print("  2-5 -> " + str(dihedral_angles["2-5"]))
    print(" 5-10 -> " + str(dihedral_angles["5-10"]))
    print("10-20 -> " + str(dihedral_angles["10-20"]))
    print("  >20 -> " + str(dihedral_angles[">20"]))


def calcJCoup(param_set, calced, experimental, Jcoup_type):
    """Calculates J-coupling values from dihedral angles
       note: all angles must be in radian"""
    JCoup_calced    = {}

    if param_set == '1':
        my_karplus = Jcoup_dict1
    elif param_set == '2':
        my_karplus = Jcoup_dict2
    elif param_set == '3':
        my_karplus = Jcoup_dict3

    A     = my_karplus['A']
    B     = my_karplus['B']
    C     = my_karplus['C']
    THETA = my_karplus['THETA']

    for record in experimental:  # resnums
        J = 0

        for my_dict in calced:  # lists (with models as dicts)
            phi = my_dict[record.resnum]

            J += (A[Jcoup_type] * (math.cos(phi + THETA[Jcoup_type])) ** 2 +
                  B[Jcoup_type] * math.cos(phi + THETA[Jcoup_type]) +
                  C[Jcoup_type])

        JCoup_calced[record.resnum] = J / len(calced)

    model_data_list = []
    model_data_dict = {}

    for Jcoup_dict in calced:   # model
        for record in experimental:
            phi = Jcoup_dict[record.resnum]

            J = (
                A[Jcoup_type] * (math.cos(phi + THETA[Jcoup_type])) ** 2 +
                B[Jcoup_type] * math.cos(phi + THETA[Jcoup_type]) +
                C[Jcoup_type])

            model_data_dict[record.resnum] = J

        model_data_list.append(model_data_dict)
        model_data_dict = {}

    return JCoup_calced, model_data_list


def calcCorrel(calced, experimental):
    """Calculates correlation between calculated and experimental data
       "calced" is a dict containing values for residues (as keys)
       "experimental" is a list containing STR record objects"""
    M = [0.0, 0.0, 0.0]
    D = [0.0, 0.0]
    match_count = 0

    for i in experimental:
        exp  = i.value
        try:
            calc = calced[i.resnum]
        except KeyError:
            continue

        M[0] += calc
        M[1] += exp
        M[2] += calc * exp

        match_count += 1

    M[0] /= match_count
    M[1] /= match_count
    M[2] /= match_count

    for i in experimental:
        exp  = i.value
        try:
            calc = calced[i.resnum]
        except KeyError:
            continue

        D[0] += (calc - M[0]) ** 2
        D[1] += (exp  - M[1]) ** 2

    D[0] /= match_count
    D[0] = math.sqrt(D[0])
    D[1] /= match_count
    D[1] = math.sqrt(D[1])

    if D[0] * D[1] == 0:
        return -2
    else:
        return (M[2] - (M[0] * M[1])) / (D[0] * D[1])


def calcQValue(calced, experimental):
    """Calculates Q-value for calculated and experimental data
       "calced" is a dict containing values for residues (as keys)
       "experimental" is a list containing STR record objects"""
    D2, E2 = 0, 0

    for i in experimental:
        exp  = i.value
        try:
            calc = calced[i.resnum]
        except KeyError:
            continue

        D2 += (calc - exp) ** 2
        E2 += exp ** 2

    Q = 100 * math.sqrt(D2) / math.sqrt(E2)
    return round(Q, 6)


def calcRMSD(calced, experimental):
    """Calculates RMSD for calculated and experimental data
       "calced" is a dict containing values for residues (as keys)
       "experimental" is a list containing STR record objects"""
    D2 = 0
    match_count = 0

    for i in experimental:
        exp  = i.value
        try:
            calc = calced[i.resnum]
        except KeyError:
            continue

        D2 += (calc - exp) ** 2

        match_count += 1

    RMSD = math.sqrt(D2 / match_count)
    return round(RMSD, 6)


def pdb2coords(PDB_file):
    """Loads PDB coordinates into a dictonary, per model"""
    model_list = PDB_model.model_list

    prev_resnum = -1
    PDB_coords = {}

    for model_num, model in enumerate(model_list):
        PDB_coords[model_num] = {}

        for atom in model:
            resnum = int(atom.getResindex()) + 1
            name   = str(atom.getName())

            if resnum == prev_resnum:
                PDB_coords[model_num][resnum][name] = Vec_3D(atom.getCoords())

            else:
                PDB_coords[model_num][resnum] = {}
                PDB_coords[model_num][resnum][name] = Vec_3D(atom.getCoords())
                prev_resnum = resnum

    return PDB_coords


def makeGraph(my_path, calced, my_experimental, graph_name):
    """X axis -> residue numbers, Y axis -> values
       "calced" is a dict containing values for residues (as keys)
       "experimental" is a list containing STR record objects"""
    experimental = copy.deepcopy(my_experimental)

    min_calc = min(calced.values())
    max_calc = max(calced.values())

    exp_values = []

    for record in experimental:
        exp_values.append(record.value)

    min_exp = min(exp_values)
    max_exp = max(exp_values)
    miny    = min(min_calc, min_exp)               # get minimum value
    maxy    = max(max_calc, max_exp)               # get maximum value

    exp_line, calc_line = [], []

    for k in range(0, max(calced.keys()) + 0):  # fetch data from arguments
        if k in list(calced.keys()):
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

    plt.figure(figsize=(10, 5), dpi=80)

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
    plt.tight_layout(pad=1.08)
    plt.savefig(my_path + "/" + graph_name, format="svg")
    plt.close()


def makeCorrelGraph(my_path, calced, experimental, graph_name):
    """X axis -> experimental values, Y axis -> calculated values
       "calced" is a dict containing values for residues (as keys)
       "experimental" is a list containing STR record objects"""
    min_calc = min(calced.values())
    max_calc = max(calced.values())

    exp_values = []
    for record in experimental:
        exp_values.append(record.value)

    min_exp = min(exp_values)
    max_exp = max(exp_values)
    miny    = min(min_calc, min_exp)             # get minimum value
    maxy    = max(max_calc, max_exp)             # get maximum value

    exp_line, calc_line = [], []

    for i, j in enumerate(calced.keys()):        # fetch data from arguments
        calc = calced[j]
        exp  = experimental[i].value

        exp_line.append(exp)
        calc_line.append(calc)

    diag = []

    for i in np.arange(miny, maxy * 1.42, 0.1):  # draw graph diagonal
        diag.append(i)

    plt.figure(figsize=(6, 5), dpi=80)
    plt.plot(diag, diag, linewidth=2.0, color='red')
    plt.plot(exp_line, calc_line, 'bo')
    plt.axis([miny, maxy, miny, maxy])
    plt.xlabel('experimental')
    plt.ylabel('calculated')
    plt.tight_layout(pad=1.08)
    plt.savefig(my_path + "/" + graph_name, format="svg")
    plt.close()


def modCorrelGraph(my_path, correl, avg_corr, model_corrs, corr_graph_name):
    """Y axis -> correlation values
       X axis -> ensemble correlation, model avg. correlation,
                 per modeel correlation
       parameter 'model_corrs' is a list containing per model correlation values
       """
    plt.figure(figsize=(6, 5), dpi=80)

    plt.plot(list(range(0, len(model_corrs))), [correl] * len(model_corrs),
             linewidth=2.0, color='green', label='Ensemble corr.')
    plt.plot(list(range(0, len(model_corrs))), [avg_corr] * len(model_corrs),
             linewidth=2.0, color='red', label='Avg. corr. per model')
    plt.plot(list(range(0, len(model_corrs))), sorted(model_corrs),
             linewidth=2.0, color='blue', label='Corr. per model')

    plt.legend(loc='lower left')
    plt.axis([-1, len(model_corrs), 0, 1])
    plt.xlabel('models (worse to best)')
    plt.ylabel('correlation')
    plt.tight_layout(pad=1.08)
    plt.savefig(my_path + "/" + corr_graph_name, format="svg")
    plt.close()


def makeNOEHist(my_path, violations):
    plt.figure(figsize=(6, 5), dpi=80)
    n_groups = len(violations)

    means_men = [
        violations['0-0.5'], violations['0.5-1'], violations['1-1.5'],
        violations['1.5-2'], violations['2-2.5'], violations['2.5-3'],
        violations['3<']
    ]

    ticks = ['0-0.5', '0.5-1', '1-1.5', '1.5-2', '2-2.5', '2.5-3', '3<']
    index = np.arange(n_groups)
    bar_width = 0.7
    plt.bar(index, means_men, bar_width, alpha=1, color='b')

    plt.xlabel("Violation (Å)")
    plt.ylabel("# of NOE distance violations")
    plt.title("NOE distance violations")
    plt.xticks(index + bar_width / 2, ticks)
    ax = plt.axes()
    ax.yaxis.grid()

    plt.tight_layout()
    plt.savefig(my_path + "/NOE_hist.svg", format="svg")
    plt.close()


def makeNMRPrideGraph(my_path, graph_data, avg_score):
    graph_data.sort()

    plt.figure(figsize=(6, 5), dpi=80)
    plt.plot(graph_data, linewidth=2.0, color='blue', label='Model scores')
    plt.plot(list(range(0, len(graph_data))), [avg_score] * len(graph_data),
             linewidth=2.0, color='green', label='Average score')
    plt.axis([-1, len(graph_data), 0, 1])
    plt.xlabel('models by score (worse to best)')
    plt.ylabel('PRIDE-NMR score')
    plt.title("PRIDE-NMR scores")
    plt.tight_layout()
    plt.legend(loc='lower left')
    ax = plt.axes()
    ax.yaxis.grid()
    plt.savefig(my_path + "/PRIDE-NMR_score.svg", format="svg")
    plt.close()
