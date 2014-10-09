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
import csx_libs.output  as csx_out

version = "1.0"


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

print("Job started with ID:", my_id)

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
pdb_models  = os.listdir("temp")                # list of models (PDB)


#-------------------------  Write file data to HTML   ------------------------#
csx_out.writeFileTable(my_path, args, my_PDB, my_id, model_count)


#------------------------  Read  and parse STR file   ------------------------#
parsed = csx_func.parseSTR(args.STR_file)


#-----------------------------  RDC calculation  -----------------------------#
# get RDC lists from STR file, each list item contains a list of record objects
RDC_lists = csx_func.get_RDC_lists(parsed.value)

for list_num, RDC_dict in enumerate(RDC_lists):

    # Pales call, results output file "pales.out"
    csx_func.callPalesOn(pdb_models, RDC_dict, args.lc_model, args.R)

    csx_out.writeRDC_table_open(my_path, "Chemical shift list", list_num + 1)

    for RDC_type in RDC_dict.keys():
        print("RDC list", list_num + 1, RDC_type)

        # get averaged RDC values -> averageRDC[residue] = value
        averageRDC = csx_func.avgPalesRDCs("pales.out", RDC_type)

        # removing records from other RDC types
        my_averageRDC = {}

        for record in RDC_dict[RDC_type]:
            my_averageRDC[record.resnum1] = averageRDC[record.resnum1]

        correl  = csx_func.calcCorrel(my_averageRDC, RDC_dict[RDC_type])
        q_value = csx_func.calcQValue(my_averageRDC, RDC_dict[RDC_type])
        rmsd    = csx_func.calcRMSD(my_averageRDC, RDC_dict[RDC_type])

        print("Correl: ", correl)
        print("Q-val:  ", q_value)
        print("RMSD:   ", rmsd)
        print()

        graph_name = str(list_num + 1) + "_RDC_" + RDC_type + ".svg"
        csx_func.makeGraph(my_path, my_averageRDC, RDC_dict[RDC_type],
                           graph_name)

        corr_graph_name = str(list_num + 1) + "_RDC_corr_" + RDC_type + ".svg"
        csx_func.makeCorrelGraph(my_path, my_averageRDC, RDC_dict[RDC_type],
                                 corr_graph_name)

        csx_out.writeRDC_data(my_path, RDC_type, len(RDC_dict[RDC_type]),
                              correl, q_value, rmsd,
                              corr_graph_name, graph_name)

    csx_out.writeRDC_table_close(my_path)


#---------------------------------  S2 calc  ---------------------------------#
# parse S2 data from STR file
S2_dict = csx_func.parseS2_STR(parsed.value)

# get averaged S2 values -> S2_calced[residue] = value
if S2_dict:

    csx_out.write_table_open(my_path, "Order parameters (S<sup>2</sup>)")

    for S2_type in S2_dict.keys():
        S2_calced = csx_func.calcS2(my_PDB, S2_dict[S2_type],
                                    args.fit, args.fit_range)

        correl  = csx_func.calcCorrel(S2_calced, S2_dict[S2_type])
        q_value = csx_func.calcQValue(S2_calced, S2_dict[S2_type])
        rmsd    = csx_func.calcRMSD(S2_calced, S2_dict[S2_type])

        print("Correl: ", correl)
        print("Q-val:  ", q_value)
        print("RMSD:   ", rmsd)
        print()

        graph_name = "RDC_" + S2_type + ".svg"
        csx_func.makeGraph(my_path, S2_calced, S2_dict[S2_type], graph_name)

        corr_graph_name = "S2_corr_" + S2_type + ".svg"
        csx_func.makeCorrelGraph(my_path, S2_calced, S2_dict[S2_type],
                                 corr_graph_name)

        csx_out.write_table_data(my_path, S2_type,
                                 len(S2_dict[S2_type]),
                                 correl, q_value, rmsd,
                                 corr_graph_name, graph_name)

    csx_out.write_table_close(my_path)


#-----------------------------  J-coupling calc  -----------------------------#
# parse J-coupling data from STR file
Jcoup_dict  = csx_func.parseJcoup_STR(parsed.value)

if Jcoup_dict:
    dihed_lists = csx_func.calcDihedAngles(my_PDB)
    csx_out.write_table_open(my_path, "Copling constants (J-coupling)")

    for Jcoup_type in Jcoup_dict.keys():

        JCoup_calced = csx_func.calcJCoup(dihed_lists, Jcoup_dict[Jcoup_type],
                                          Jcoup_type)

        correl  = csx_func.calcCorrel(JCoup_calced, Jcoup_dict[Jcoup_type])
        q_value = csx_func.calcQValue(JCoup_calced, Jcoup_dict[Jcoup_type])
        rmsd    = csx_func.calcRMSD(JCoup_calced, Jcoup_dict[Jcoup_type])

        print("Correl: ", correl)
        print("Q-val:  ", q_value)
        print("RMSD:   ", rmsd)
        print()

        graph_name = "JCoup_" + Jcoup_type + ".svg"
        csx_func.makeGraph(my_path, JCoup_calced, Jcoup_dict[Jcoup_type],
                           graph_name)

        corr_graph_name = "JCoup_corr_" + Jcoup_type + ".svg"
        csx_func.makeCorrelGraph(my_path, JCoup_calced, Jcoup_dict[Jcoup_type],
                                 corr_graph_name)

        csx_out.write_table_data(my_path, Jcoup_type,
                                 len(Jcoup_dict[Jcoup_type]),
                                 correl, q_value, rmsd,
                                 corr_graph_name, graph_name)

    csx_out.write_table_close(my_path)


#----------------------------  Chemical shift calc  --------------------------#
# parse J-coupling data from STR file
ChemShift_lists = csx_func.parseChemShift_STR(parsed.value)

if ChemShift_lists:
    CS_calced = csx_func.callShiftxOn(pdb_models)
    csx_out.writeRDC_table_open(my_path, "Chemical shift list", list_num + 1)

    for list_num, CS_list in enumerate(ChemShift_lists):

        for CS_type in CS_list.keys():
            exp_dict = {}

            for record in CS_list[CS_type]:
                exp_dict[record.resnum] = CS_calced[CS_type][record.resnum]

            correl  = csx_func.calcCorrel(exp_dict, CS_list[CS_type])
            q_value = csx_func.calcQValue(exp_dict, CS_list[CS_type])
            rmsd    = csx_func.calcRMSD(exp_dict, CS_list[CS_type])

            print("CHEM SHIFT", CS_type)
            print("Correl: ", correl)
            print("Q-val:  ", q_value)
            print("RMSD:   ", rmsd)
            print()

            graph_name = str(list_num + 1) + "_CS_" + CS_type + ".svg"
            csx_func.makeGraph(my_path, exp_dict, CS_list[CS_type],
                               graph_name)

            corr_graph_name = str(list_num + 1) + "_CS_corr_" + CS_type + ".svg"
            csx_func.makeCorrelGraph(my_path, exp_dict, CS_list[CS_type],
                                     corr_graph_name)

            csx_out.writeRDC_data(my_path, CS_type, len(CS_list[CS_type]),
                                  correl, q_value, rmsd,
                                  corr_graph_name, graph_name)

    csx_out.writeRDC_table_close(my_path)




# if args.PDB_fetch:
#     os.remove(my_PDB)


