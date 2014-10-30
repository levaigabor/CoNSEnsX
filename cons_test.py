#!/usr/bin/python
# -*- coding: utf-8 -*-

# standard modules
from __future__ import print_function
from multiprocessing import Process, Pipe
import os
import string
import random
import time
import subprocess
import ast

import math

# own modules
import csx_libs.calc    as csx_calc
import csx_libs.misc    as csx_misc
import csx_libs.methods as csx_func
import csx_libs.objects as csx_obj
import csx_libs.output  as csx_out

version = "1.0"

ts = time.time()

#------------------  Setting up parser and parse CLI arguments  --------------#
parser = csx_misc.createParser()            # get parser from module
args   = parser.parse_args()                # parsing CLI arguments

if args.c:                                  # show credit
    print(csx_misc.cred)
    raise SystemExit

if args.help:                               # show help
    print(csx_misc.help_text)
    raise SystemExit

if args.STR_file and (args.PDB_file or args.PDB_fetch):
    pass
else:
    print("missing file")
    parser.print_usage()
    raise SystemExit

if args.PDB_fetch and args.PDB_file:
    print("Double input for PDB data, deselect one of the input sources")
    raise SystemExit

my_PDB = csx_func.get_PDB(args)


#------------------------  Setting up working directory  ---------------------#
def getID(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

my_id = getID()
my_path = "calculations/" + my_id + '/'

print("Job started with ID: \033[0;35m" + my_id + "\033[0m")

if not os.path.exists(my_path):                # create working folder
    os.makedirs(my_path)
else:
    for f in os.listdir(my_path):              # clean working folder
        os.remove(my_path + '/' + f)


#----------------  Setting up output files & write header data  --------------#
# csx_out.writeHeaderTXT(my_path, args, version)

csx_out.writeHeaderHTML(my_path, version)


#-------------------------  Making temporary folder   ------------------------#
if not os.path.exists("temp"):                 # create temp folder
    os.makedirs("temp")
else:
    for f in os.listdir("temp"):               # clean temp folder
        os.remove("temp/" + f)

csx_func.pdb_cleaner(my_PDB)                   # bringing PDB to format
model_count = csx_func.pdb_splitter(my_PDB)    # splitting of PDB file
pdb_models  = os.listdir("temp")               # list of models (PDB)


#-------------------------  Write file data to HTML   ------------------------#
csx_out.writeFileTable(my_path, args, my_PDB, my_id, model_count)


#------------------------  Read  and parse STR file   ------------------------#
parsed = csx_func.parseSTR(args.STR_file)

def f(conn, restain_file):
    print("PARSER PROCESS STARTED")
    # restraints_parsed = csx_func.parseSTR(args.XPLOR_file)


    restraints_parsed = subprocess.check_output(["pypy", "pypyParse.py",
                                                args.XPLOR_file])

    restraints_parsed = ast.literal_eval(restraints_parsed)

    conn.send(restraints_parsed)

    print("PARSER PROCESS READY")
    conn.close()


if args.XPLOR_file:
    parent_conn, child_conn = Pipe()
    p = Process(target=f, args=(child_conn, args.XPLOR_file,))
    p.start()


#-----------------------------  RDC calculation  -----------------------------#
# RDC_lists = csx_func.get_RDC_lists(parsed.value)

# if RDC_lists:
#     csx_calc.calcRDC(RDC_lists, pdb_models, my_path, args)


#---------------------------------  S2 calc  ---------------------------------#
# S2_dict = csx_func.parseS2_STR(parsed.value)

# if S2_dict:
#     csx_calc.calcS2(S2_dict, my_PDB, my_path, args)


#-----------------------------  J-coupling calc  -----------------------------#
# Jcoup_dict  = csx_func.parseJcoup_STR(parsed.value)

# if Jcoup_dict:
#     csx_calc.calcJCouplings(Jcoup_dict, my_PDB, my_path)


#----------------------------  Chemical shift calc  --------------------------#
# ChemShift_lists = csx_func.parseChemShift_STR(parsed.value)

# if ChemShift_lists:
#     csx_calc.calcChemShifts(ChemShift_lists, pdb_models, my_path)




if args.XPLOR_file:
    saveShifts = parent_conn.recv()
    p.join()

    restraints = []

    # parse data to restraint objects returned from pypy process
    for data in saveShifts:
        restraints.append(csx_obj.Restraint_Record(
                          data[0], data[1], data[2], data[3], data[4], data[5]))

    PDB_coords    = csx_func.parse2dicts(args.PDB_file)
    prev_id       = -1
    avg_distances = {}
    str_distaces  = {}

    for restraint in restraints:

        ######################### TEMPONARY WORKAROUND #########################
        if restraint.atom_ID1 == "ME" or restraint.atom_ID2 == "ME":
            continue
        if restraint.atom_ID1 == "MG" or restraint.atom_ID2 == "MG":
            continue
        if restraint.atom_ID1 == "MB" or restraint.atom_ID2 == "MB":
            continue
        if restraint.atom_ID1 == "MD" or restraint.atom_ID2 == "MD":
            continue
        if restraint.atom_ID1 == "MD1" or restraint.atom_ID2 == "MD1":
            continue
        if restraint.atom_ID1 == "MD2" or restraint.atom_ID2 == "MD2":
            continue
        if restraint.atom_ID1 == "MG1" or restraint.atom_ID2 == "MG1":
            continue
        if restraint.atom_ID1 == "MG2" or restraint.atom_ID2 == "MG2":
            continue
        ######################### TEMPONARY WORKAROUND #########################

        curr_id = int(restraint.curr_distID)

        if prev_id == curr_id:
            model_avg_dist = csx_func.getModelAvgDistance(PDB_coords,
                                                          restraint.seq_ID1,
                                                          restraint.atom_ID1,
                                                          restraint.seq_ID2,
                                                          restraint.atom_ID2)

            avg_distances[curr_id].append(model_avg_dist)

        else:
            prev_id = curr_id
            avg_distances[curr_id] = []
            str_distaces[curr_id] = restraint.dist_max

            model_avg_dist = csx_func.getModelAvgDistance(PDB_coords,
                                                          restraint.seq_ID1,
                                                          restraint.atom_ID1,
                                                          restraint.seq_ID2,
                                                          restraint.atom_ID2)

            avg_distances[curr_id].append(model_avg_dist)


    # averaging over the same restraint ID data
    for key in avg_distances.keys():
        avg = 0.0

        for distance in avg_distances[key]:
            avg += math.pow(float(distance), -6)

        avg_distances[key] = math.pow(avg / len(avg_distances[key]), -1.0/6)


    avg_dist_keys = avg_distances.keys()
    avg_dist_keys.sort()
    violations = 0

    for key in avg_dist_keys:
        if avg_distances[key] > str_distaces[key]:
            # print(key, "NOE VIOLATION at")
            violations += 1


    print("Total # of violations:", violations)



te = time.time()
print("total runtime", te-ts)
print("Your ID was: \033[0;35m" + my_id + "\033[0m")
