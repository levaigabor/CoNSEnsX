#!/usr/bin/python
# -*- coding: utf-8 -*-

# standard modules
from __future__ import print_function
import os

# installed modules
import nmrpystar

# own modules
import csx_libs.misc    as csx_misc
import csx_libs.methods as csx_func
import csx_libs.objects as csx_obj
import csx_libs.ouput   as csx_out

version = "1.0"

parser = csx_misc.createParser()            # get parser from module
args = parser.parse_args()                  # parsing CLI arguments

if args.c:                                  # show credit
    print(csx_misc.cred)
    raise SystemExit

if args.help:                               # show help
    print(csx_misc.help_text)
    raise SystemExit

if not args.STR_file or not args.PDB_file:  # checking for input files
    print("missing file")
    parser.print_usage()
    raise SystemExit


#----------------  Setting up output files & write header data  --------------#
if "txt" in args.output_format:
    csx_out.writeHeaderTXT(args, version)

if "html" in args.output_format:
    csx_out.writeHeaderHTML(args, version)


#-------------------------  Making temporary folder   ------------------------#
if not os.path.exists("temp"):          # create temp folder
    os.makedirs("temp")
else:
    for f in os.listdir("temp"):        # clean temp folder
        os.remove("temp/" + f)

csx_func.pdb_cleaner(args.PDB_file)     # bringing PDB to format
csx_func.pdb_splitter(args.PDB_file)    # splitting of PDB file


#------------------------  Read  and parse STR file   ------------------------#
star_file = open(args.STR_file)         # open STR file
myString = ""

for line in star_file:                  # rean STR file into a string
    myString += line

parsed = nmrpystar.parse(myString)      # parsing, access data -> parsed.value

if parsed.status != 'success':          # check if parsing was successful
    print('Error during STR parsing: ', parsed)
    raise SystemExit


#-----------------------------  RDC calculation  -----------------------------#
# get RDC lists from STR file, each list item contains a list of record objects
RDC_lists  = csx_func.get_RDC_lists(parsed.value)
pdb_models = os.listdir("temp")         # list of models (PDB)

for RDC_list in RDC_lists:
    #### for RDC_type in RDC_list.keys()

    # Pales call, results output file "pales.out"
    csx_func.callPalesOn(pdb_models, RDC_list, args.lc_model, args.R)

    # get averaged RDC values -> averageRDC[residue] = value
    averageRDC = csx_func.avgPalesRDCs("pales.out")

    print("Correl: ", csx_func.calcCorrel(averageRDC, RDC_list))
    print("Q-val:  ", csx_func.calcQValue(averageRDC, RDC_list))
    print("RMSD:   ", csx_func.calcRMSD(averageRDC, RDC_list))
    csx_func.makeGraph(averageRDC, RDC_list)
    csx_func.makeCorrelGraph(averageRDC, RDC_list)

    #### for RDC_type in RDC_list.keys()


#------------------------ parse S2 data from STR file   ----------------------#
# get S2 values from STR file, each list item contains a list of record objects
try:
    saveShifts = parsed.value.saves["order_param"]
except KeyError:
    print("No S2 parameter list found")

loopShifts = saveShifts.loops[-1]
S2_records = []

for ix in range(len(loopShifts.rows)):   # fetch values from STR file
            row = loopShifts.getRowAsDict(ix)

            S2_records.append(csx_obj.S2_Record(row["Residue_seq_code"],
                                                row["Atom_name"],
                                                row["S2_value"]))


#---------------------------------  S2 calc  ---------------------------------#
# get averaged S2 values -> S2_calced[residue] = value
S2_calced = csx_func.calcS2(args.PDB_file, S2_records)

print("S2_corr:", csx_func.calcCorrel(S2_calced, S2_records))
print("S2Q-val:", csx_func.calcQValue(S2_calced, S2_records))
print("RMSD:   ", csx_func.calcRMSD(S2_calced, S2_records))
csx_func.makeGraph(S2_calced, S2_records)
csx_func.makeCorrelGraph(S2_calced, S2_records)
