#!/usr/bin/python
# -*- coding: utf-8 -*-

# standard modules
from __future__ import print_function
import os
import string
import random

# own modules
import csx_libs.misc    as csx_misc
import csx_libs.methods as csx_func
import csx_libs.objects as csx_obj
import csx_libs.output   as csx_out

version = "1.0"

#------------------------  Setting up working directory  ---------------------#
def getID(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

my_id = getID()
my_path = "calculations/" + my_id + '/'

if not os.path.exists(my_path):          # create temp folder
    os.makedirs(my_path)
else:
    for f in os.listdir(my_path):        # clean temp folder
        os.remove(my_path + '/' + f)


#------------------  Setting up parser and parse CLI arguments  --------------#
parser = csx_misc.createParser()            # get parser from module
args   = parser.parse_args()                # parsing CLI arguments

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
# if "txt" in args.output_format:
#     csx_out.writeHeaderTXT(my_path, args, version)

if "html" in args.output_format:
    csx_out.writeHeaderHTML(my_path, version)


#-------------------------  Making temporary folder   ------------------------#
if not os.path.exists("temp"):                        # create temp folder
    os.makedirs("temp")
else:
    for f in os.listdir("temp"):                      # clean temp folder
        os.remove("temp/" + f)

csx_func.pdb_cleaner(args.PDB_file)                   # bringing PDB to format
model_count = csx_func.pdb_splitter(args.PDB_file)    # splitting of PDB file


#-------------------------  Write file data to HTML   ------------------------#
csx_out.writeFileTable(my_path, args, my_id, model_count)


#------------------------  Read  and parse STR file   ------------------------#
parsed = csx_func.parseSTR(args.STR_file)


#-----------------------------  RDC calculation  -----------------------------#
# get RDC lists from STR file, each list item contains a list of record objects
RDC_lists  = csx_func.get_RDC_lists(parsed.value)
pdb_models = os.listdir("temp")         # list of models (PDB)

prev_RDC_dict = {}

for list_num, RDC_dict in enumerate(RDC_lists):

    # Pales call, results output file "pales.out"
    csx_func.callPalesOn(pdb_models, RDC_dict, args.lc_model, args.R)

    csx_out.writeRDC_table_open(my_path, list_num + 1)

    for RDC_type in RDC_dict.keys():
        print("RDC list", list_num + 1, RDC_type)

        # get averaged RDC values -> averageRDC[residue] = value
        averageRDC = csx_func.avgPalesRDCs("pales.out", RDC_type)

        # removing records from other RDC types
        my_averageRDC = {}

        for record in RDC_dict[RDC_type]:
            my_averageRDC[record.resnum1] = averageRDC[record.resnum1]

        print("Correl: ", csx_func.calcCorrel(my_averageRDC, RDC_dict[RDC_type]))
        print("Q-val:  ", csx_func.calcQValue(my_averageRDC, RDC_dict[RDC_type]))
        print("RMSD:   ", csx_func.calcRMSD(my_averageRDC, RDC_dict[RDC_type]))
        print()
        csx_func.makeGraph(my_averageRDC, RDC_dict[RDC_type],
                           "RDC_" + RDC_type)
        csx_func.makeCorrelGraph(my_averageRDC, RDC_dict[RDC_type],
                                 "RDC_correlation_" + RDC_type)

    csx_out.writeRDC_table_close(my_path)

#---------------------------------  S2 calc  ---------------------------------#
# parse S2 data from STR file
S2_dict = csx_func.parseS2_STR(parsed.value)

# get averaged S2 values -> S2_calced[residue] = value
if S2_dict:
    for S2_type in S2_dict.keys():
        S2_calced = csx_func.calcS2(args.PDB_file, S2_dict[S2_type],
                                    args.fit, args.fit_range)

        print("S2_corr:", csx_func.calcCorrel(S2_calced, S2_dict[S2_type]))
        print("S2Q-val:", csx_func.calcQValue(S2_calced, S2_dict[S2_type]))
        print("RMSD:   ", csx_func.calcRMSD(S2_calced, S2_dict[S2_type]))
        print()
        # csx_func.makeGraph(S2_calced, S2_dict[S2_type])
        # csx_func.makeCorrelGraph(S2_calced, S2_dict[S2_type])


#-----------------------------  J-coupling calc  -----------------------------#
# parse J-coupling data from STR file
Jcoup_dict = csx_func.parseJcoup_STR(parsed.value)
avg_dict   = csx_func.calcDihedAngles(args.PDB_file)

if Jcoup_dict:
    for Jcoup_type in Jcoup_dict.keys():

        JCoup_calced = csx_func.calcJCoup(avg_dict, Jcoup_dict[Jcoup_type],
                                          Jcoup_type)

        print(Jcoup_type + "_corr:", csx_func.calcCorrel(JCoup_calced,
                                                    Jcoup_dict[Jcoup_type]))
        print("Q-val:  ", csx_func.calcQValue(JCoup_calced, Jcoup_dict[Jcoup_type]))
        print("RMSD:   ", csx_func.calcRMSD(JCoup_calced, Jcoup_dict[Jcoup_type]))
        # csx_func.makeGraph(JCoup_calced, Jcoup_dict[Jcoup_type],
        #                    "J-coupling " + Jcoup_type)
        # csx_func.makeCorrelGraph(JCoup_calced, Jcoup_dict[Jcoup_type])











