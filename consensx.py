#!/usr/bin/python
# -*- coding: utf-8 -*-

# standard modules
from __future__ import print_function
from multiprocessing import Process, Pipe
import os
import time
import subprocess
import ast

# own modules
import csx_libs.calc    as csx_calc
import csx_libs.misc    as csx_misc
import csx_libs.methods as csx_func
import csx_libs.output  as csx_out
import csx_libs.objects as csx_obj

version  = "0.9_python"
pales    = "/home/daniel/Programme/linux/pales"
shiftx   = "/home/daniel/Programme/shiftx/shiftx"
pridedb  = "/home/daniel/Programme/prideNMR/pdb2hhbindbM"
pridenmr = "/home/daniel/Programme/prideNMR/mrhisthhbindbM"


ts = time.time()

#------------------  Setting up parser and parse CLI arguments  ---------------#
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


#----------------  Setting up working directory and results HTML  -------------#
my_id, my_path = csx_func.getID()
os.makedirs(my_path)

csx_out.writeHeaderHTML(my_path, version)


#-----------------------------  Prepare PDB data  -----------------------------#
my_PDB = csx_func.get_PDB(args)
csx_func.pdb_cleaner(my_PDB)                   # bringing PDB to format
csx_func.get_model_list(my_PDB)
model_count = csx_func.pdb_splitter(my_path, my_PDB)

pdb_models = []                                # list of models (PDB)
for file in os.listdir(my_path):
    if file.endswith(".pdb"):
        pdb_models.append(file)

csx_func.calcPeptideBonds(my_PDB)

#------------------------  Read  and parse STR file   -------------------------#
parsed = csx_func.parseSTR(args.STR_file)

def pypyProcess(conn, restain_file):
    """function called in seperate process to parse mr file"""
    pypy_command      = ["pypy", "pypyParse.py", args.NOE_file]
    restraints_parsed = subprocess.check_output(pypy_command)
    restraints_parsed = ast.literal_eval(restraints_parsed)
    conn.send(restraints_parsed)
    conn.close()

if args.NOE_file:
    parent_conn, child_conn = Pipe()
    p = Process(target = pypyProcess, args = (child_conn, args.NOE_file,))
    p.start()
#------------------------  NOE distance violation calc  -----------------------#
# SEPARATE HERE TO JOIN PROCESS LATER xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# if args.NOE_file:
    saveShifts = parent_conn.recv()
    p.join()
    NOE_violations = csx_calc.calcNOEviolations(args, saveShifts,
                                                my_path, args.r3_averaging)
    PRIDE_data = csx_calc.calcNMR_Pride(pdb_models, my_path)
    csx_out.writeFileTable(my_path, args, my_PDB, my_id,
                           model_count, args.NOE_file,
                           csx_obj.Restraint_Record.getRestraintCount())
    csx_out.write_bottom_table(my_path, NOE_violations, PRIDE_data)
# SEPARATE HERE TO JOIN PROCESS LATER xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
else:
    csx_out.writeFileTable(my_path, args, my_PDB, my_id, model_count)


#-----------------------------  RDC calculation  ------------------------------#
RDC_lists = csx_func.get_RDC_lists(parsed.value)

if RDC_lists:
    csx_calc.calcRDC(RDC_lists, pdb_models, my_path, args)


#---------------------------------  S2 calc  ----------------------------------#
S2_dict = csx_func.parseS2_STR(parsed.value)

if S2_dict:
    csx_calc.calcS2(S2_dict, my_PDB, my_path, args)


#-----------------------------  J-coupling calc  ------------------------------#
Jcoup_dict  = csx_func.parseJcoup_STR(parsed.value)

if Jcoup_dict:
    csx_calc.calcJCouplings(Jcoup_dict, my_PDB, my_path)


#----------------------------  Chemical shift calc  ---------------------------#
ChemShift_lists = csx_func.parseChemShift_STR(parsed.value)

if ChemShift_lists:
    csx_calc.calcChemShifts(ChemShift_lists, pdb_models, my_path)


te = time.time()
print("total runtime", te-ts)
print("Your ID was: \033[0;35m" + my_id + "\033[0m")
