#!/usr/bin/env python3

"""
              _____      _   _  _____ ______          __   __
             / ____|    | \ | |/ ____|  ____|         \ \ / /
            | |     ___ |  \| | (___ | |__   _ __  ___ \ V /
            | |    / _ \| . ` |\___ \|  __| | '_ \/ __| > <
            | |___| (_) | |\  |____) | |____| | | \__ \/ . \
             \_____\___/|_| \_|_____/|______|_| |_|___/_/ \_\

     Compliance of NMR-derived Structural Ensembles with experimental data

Authors:    Zoltán Gáspári, Dániel Dudola
Fork me at: https://github.com/derPuntigamer/CoNSEnsX
"""

# standard modules
import os
import time

# own modules
import csx_libs.calc    as csx_calc
import csx_libs.parser  as csx_parser
import csx_libs.methods as csx_func
import csx_libs.output  as csx_out
import csx_libs.objects as csx_obj

import csx_sel.selection as csx_sel

version = 0.4

ts = time.time()

abspath = os.path.abspath(__file__)
dirname = os.path.dirname(abspath)

csx_func.check_3rd_party(dirname)

#------------------  Setting up parser and parse CLI arguments  ---------------#
parser = csx_parser.createParser()            # get parser from module
args   = parser.parse_args()                # parsing CLI arguments

if args.c:                                  # show credit
    print(csx_parser.cred)
    raise SystemExit

if args.help:                               # show help
    print(csx_parser.help_text)
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
my_id, my_path = csx_func.getID(args)
os.makedirs(my_path)
csx_obj.CSV_buffer.working_dir = my_path

csx_out.writeHeaderHTML(my_path, version)


#-----------------------------  Prepare PDB data  -----------------------------#
my_PDB = csx_func.get_PDB(args)
#csx_func.pdb_cleaner(my_PDB)                   # bringing PDB to format
csx_func.get_model_list(my_PDB, args.PDB_file)
model_count = csx_func.pdb_splitter(my_path, my_PDB)

pdb_models = []                                # list of models (PDB)
for file in os.listdir(my_path):
    if file.endswith(".pdb"):
        pdb_models.append(file)

pdb_models = csx_func.natural_sort(pdb_models)

if args.a:
    csx_func.calcPeptideBonds(my_PDB)
    csx_func.calcNH_Angles(my_PDB)


#------------------------  Read  and parse STR file   -------------------------#
parsed = csx_func.parseSTR(args.STR_file)

if args.NOE_file:
    saveShifts = csx_func.getNOE(args.NOE_file)
    NOE_violations = csx_calc.calcNOEviolations(args, saveShifts,
                                                my_path, args.r3_averaging)
    PRIDE_data = csx_calc.calcNMR_Pride(pdb_models, my_path)
    csx_out.writeFileTable(my_path, args, my_PDB, my_id,
                           model_count, args.NOE_file,
                           csx_obj.Restraint_Record.getRestraintCount())
    csx_out.write_bottom_table(my_path, NOE_violations, PRIDE_data)
else:
    csx_out.writeFileTable(my_path, args, my_PDB, my_id, model_count)


#-----------------------------  RDC calculation  ------------------------------#
RDC_lists = csx_func.get_RDC_lists(parsed.value)

if RDC_lists:
    csx_calc.calcRDC(RDC_lists, pdb_models, my_path, args)


#---------------------------------  S2 calc  ----------------------------------#
S2_dict = csx_func.parseS2_STR(parsed.value)

if S2_dict:
    csx_calc.calcS2(S2_dict, my_path, args)

S2_sidechain = csx_func.parse_sidechain_S2_STR(parsed.value)

if S2_sidechain:
    csx_calc.calcS2_sidechain(S2_sidechain, my_path, args)


#-----------------------------  J-coupling calc  ------------------------------#
Jcoup_dict = csx_func.parseJcoup_STR(parsed.value)

if Jcoup_dict:
    csx_calc.calcJCouplings(args.d, Jcoup_dict, my_PDB, my_path)


#----------------------------  Chemical shift calc  ---------------------------#
ChemShift_lists = csx_func.parseChemShift_STR(parsed.value)

if ChemShift_lists:
    csx_calc.calcChemShifts(ChemShift_lists, pdb_models, my_path)

csx_obj.CSV_buffer.writeCSV()
csx_out.close_HTML(my_path)

te = time.time()
print("total runtime", te-ts)
print("Your ID was: \033[0;35m" + my_id + "\033[0m")


# user_sel = [
#     #["RDC", 1, "0_C_CA", 0.5],
#     ["RDC", 1, "0_N_H", 0.5],
#     #["S2", "N", 0.25],
#     ["JCoup", "3JHNHA", 0.5],
#     ["ChemShift", "CA", 0.9]
# ]

# csx_sel.selection_on("RDC", "correlation", pdb_models, RDC_lists, user_sel,
#                      args=args, S2_dict=S2_dict, Jcoup_dict=Jcoup_dict,
#                      ChemShifts=ChemShift_lists,
#                      S2_type="N", overdrive=2)

