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

csx_func.check_3rd_party()

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
csx_func.get_model_list(my_PDB)
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
Jcoup_dict  = csx_func.parseJcoup_STR(parsed.value)

if Jcoup_dict:
    csx_calc.calcJCouplings(args.d, Jcoup_dict, my_PDB, my_path)


#----------------------------  Chemical shift calc  ---------------------------#
ChemShift_lists = csx_func.parseChemShift_STR(parsed.value)

if ChemShift_lists:
    csx_calc.calcChemShifts(ChemShift_lists, pdb_models, my_path)

csx_obj.CSV_buffer.writeCSV()

te = time.time()
print("total runtime", te-ts)
print("Your ID was: \033[0;35m" + my_id + "\033[0m")




user_sel = [["RDC", 2, "0_N_H", 1], ["RDC", 1, "0_N_H", 0.5], ["RDC", 3, "0_N_H", 0.5]]


#print(RDC_lists)

best_num     = csx_sel.RDC_modell_data.get_best_model(user_sel, RDC_lists)
in_selection = [best_num] # int number

first_run = True
prev_best = -2



while True:
    model_scores = {}

    for num, pdb in enumerate(pdb_models):
        if num in in_selection:
            continue

        divide_by = 0.0
        pdb_sel = [num]

        for selected in in_selection:
            pdb_sel.append(selected)

        print("current selection: ", pdb_sel)

        for sel_data in user_sel:
            if sel_data[0] != "RDC":
                continue

            RDC_num    = sel_data[1]
            RDC_type   = sel_data[2]
            RDC_weight = sel_data[3]

            #my_data = csx_sel.RDC_modell_data.RDC_data[RDC_num][RDC_type][num]
            #pdb_sel.append(my_data)

            print(pdb_sel)

            averageRDC = csx_sel.averageRDCs_on(pdb_sel, RDC_num, RDC_type)
            correl  = csx_func.calcCorrel(averageRDC, RDC_lists[RDC_num - 1][RDC_type])

            if num in model_scores.keys():
                model_scores[num] += correl * RDC_weight
            else:
                model_scores[num] = correl * RDC_weight

            divide_by += RDC_weight

    best_num = -1
    best_val = -1

    for num in model_scores.keys():
        if model_scores[num] / divide_by > best_val:
            best_val = model_scores[num] / divide_by
            best_num = num

    print("prev best: " + str(prev_best) + ", current best: " + str(best_val))

    if best_val > prev_best:
        prev_best = best_val
        in_selection.append(best_num)
    else:
        #in_selection = [x+1 for x in in_selection]
        in_selection.sort()
        print("numbered as in PDB file:\n", in_selection)
        raise SystemExit
