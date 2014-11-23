#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function
import __main__
import math
import subprocess
import os

import methods as csx_func
import objects as csx_obj
import output  as csx_out

import matplotlib.pyplot as plt

def calcRDC(RDC_lists, pdb_models, my_path, args):
    for list_num, RDC_dict in enumerate(RDC_lists):

        # Pales call, results output file "pales.out"
        csx_func.callPalesOn(my_path, pdb_models, RDC_dict, args.lc_model, args.R)

        csx_out.writeRDC_table_open(my_path, "RDC list", list_num + 1)

        for RDC_type in RDC_dict.keys():
            print("RDC list", list_num + 1, RDC_type)

            # get averaged RDC values -> averageRDC[residue] = value
            pales_out = my_path + "pales.out"
            averageRDC, model_data = csx_func.avgPalesRDCs(pales_out, RDC_type)

            model_corrs = []

            for model in model_data:
                model_corrs.append(csx_func.calcCorrel(model, RDC_dict[RDC_type]))

            avg_model_corr = sum(model_corrs) / len(model_corrs)

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

            mod_corr_graph_name = (str(list_num + 1) + "_RDC_mod_corr_" +
                                   RDC_type + ".svg")
            csx_func.modCorrelGraph(my_path, correl, avg_model_corr, model_corrs,
                                    mod_corr_graph_name)

            csx_out.writeRDC_data(my_path, RDC_type, len(RDC_dict[RDC_type]),
                                  correl, q_value, rmsd,
                                  corr_graph_name, graph_name, mod_corr_graph_name)

        csx_out.writeRDC_table_close(my_path)


def calcS2(S2_dict, my_PDB, my_path, args):
    csx_out.write_table_open(my_path, "Order parameters (S<sup>2</sup>)")

    for S2_type in S2_dict.keys():
        S2_calced = csx_func.calcS2(my_PDB, S2_dict[S2_type],
                                    args.fit, args.fit_range)

        correl  = csx_func.calcCorrel(S2_calced, S2_dict[S2_type])
        q_value = csx_func.calcQValue(S2_calced, S2_dict[S2_type])
        rmsd    = csx_func.calcRMSD(S2_calced, S2_dict[S2_type])

        print("N-H Order Parameters")
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


def calcJCouplings(Jcoup_dict, my_PDB, my_path):
    dihed_lists = csx_func.calcDihedAngles(my_PDB)
    csx_out.write_table_open(my_path, "Coupling constants (J-coupling)")

    for Jcoup_type in Jcoup_dict.keys():

        JCoup_calced, model_data = csx_func.calcJCoup(dihed_lists,
                                                      Jcoup_dict[Jcoup_type],
                                                      Jcoup_type)

        model_corrs = []

        for model in model_data:
            model_corrs.append(csx_func.calcCorrel(model,
                                                   Jcoup_dict[Jcoup_type]))

        avg_model_corr = sum(model_corrs) / len(model_corrs)

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

        mod_corr_graph_name = "JCoup_mod_corr_" + Jcoup_type + ".svg"
        csx_func.modCorrelGraph(my_path, correl, avg_model_corr, model_corrs,
                                mod_corr_graph_name)

        csx_out.write_table_data(my_path, Jcoup_type,
                                 len(Jcoup_dict[Jcoup_type]),
                                 correl, q_value, rmsd,
                                 corr_graph_name, graph_name, mod_corr_graph_name)

    csx_out.write_table_close(my_path)


def calcChemShifts(ChemShift_lists, pdb_models, my_path):
    CS_calced, model_data = csx_func.callShiftxOn(my_path, pdb_models)

    for list_num, CS_list in enumerate(ChemShift_lists):

        csx_out.writeRDC_table_open(my_path, "Chemical shift list",
                                    list_num + 1)

        for CS_type in CS_list.keys():
            model_corrs = []

            for model in model_data:
                inner_exp_dict = {}

                for record in CS_list[CS_type]:
                    inner_exp_dict[record.resnum] = model[CS_type][record.resnum]

                model_corrs.append(csx_func.calcCorrel(inner_exp_dict,
                                                       CS_list[CS_type]))

            avg_model_corr = sum(model_corrs) / len(model_corrs)

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

            mod_corr_graph_name = "CS_mod_corr_" + CS_type + ".svg"
            csx_func.modCorrelGraph(my_path, correl, avg_model_corr, model_corrs,
                                mod_corr_graph_name)

            csx_out.writeRDC_data(my_path, CS_type, len(CS_list[CS_type]),
                                  correl, q_value, rmsd,
                                  corr_graph_name, graph_name,
                                  mod_corr_graph_name)

    csx_out.writeRDC_table_close(my_path)


def calcNOEviolations(args, saveShifts, my_path, r3_averaging):
    # parse data to restraint objects returned from pypy process
    for data in saveShifts:
        csx_obj.Restraint_Record(data[0], data[1], data[2], data[3],
                                 data[4], data[5], data[6], data[7])

    # fetch all restraint from class
    restraints = csx_obj.Restraint_Record.all_restraints

    PDB_coords    = csx_func.pdb2coords(args.PDB_file)
    prev_id       = -1
    avg_distances = {}
    measured_avg  = {}
    str_distaces  = {}

    for model in PDB_coords.keys():
        avg_distances[model] = {}

        for restraint in restraints:
            curr_id = int(restraint.curr_distID)
            resnum1 = restraint.seq_ID1
            atom1   = restraint.atom_ID1
            resnum2 = restraint.seq_ID2
            atom2   = restraint.atom_ID2

            # print(curr_id, resnum1, atom1, resnum2, atom2)

            atom_coord1 = PDB_coords[model][resnum1][atom1]
            atom_coord2 = PDB_coords[model][resnum2][atom2]

            distance = (atom_coord1 - atom_coord2).magnitude()

            if prev_id == curr_id:
                avg_distances[model][curr_id].append(distance)

            else:
                prev_id = curr_id
                avg_distances[model][curr_id] = []
                str_distaces[curr_id] = restraint.dist_max

                avg_distances[model][curr_id].append(distance)

    # at this point avg_distances[model][curr_id] contains distances for one
    # model and one restraint GROUP identified with "curr_distID" number

    prev_id  = -1

    for model in PDB_coords.keys():
        for restraint in restraints:
            curr_id = int(restraint.curr_distID)

            if prev_id == curr_id:
                continue
            else:
                prev_id = curr_id

            avg = 0.0

            for distance in avg_distances[model][curr_id]:
                if r3_averaging:
                    avg += math.pow(float(distance), -3)
                else:
                    avg += math.pow(float(distance), -6)

            avg /= len(avg_distances[model][curr_id])
            if r3_averaging:
                avg_distances[model][curr_id] = math.pow(avg, -1.0/3)
            else:
                avg_distances[model][curr_id] = math.pow(avg, -1.0/6)

    # at this point avg_distances[model][curr_id] contain a single (r-6)
    # averaged distance for one model and one restraint GROUP identified with
    # "curr_distID" number. Averaging is done on "in GROUP" distances

    for restraint in restraints:
        curr_id = int(restraint.curr_distID)
        avg     = 0.0

        for model in PDB_coords.keys():
            avg += math.pow(avg_distances[model][curr_id], -6)

        avg /= len(PDB_coords.keys())
        measured_avg[curr_id] = math.pow(avg, -1.0/6)

    # at this point measured_avg[curr_id] is a simple dictonary containing the
    # model averaged distances for the given "curr_distID" number

    avg_dist_keys = measured_avg.keys()
    avg_dist_keys.sort()
    violations = {"0-0.5" : 0, "0.5-1" : 0, "1-1.5" : 0,
                  "1.5-2" : 0, "2-2.5" : 0, "2.5-3" : 0, "3<" : 0}
    viol_count = 0

    diff_sum = 0

    for key in avg_dist_keys:
        if measured_avg[key] > str_distaces[key]:
            viol_count += 1
            diff = measured_avg[key] - str_distaces[key]
            diff_sum += diff

            if diff <= 0.5:
                violations["0-0.5"] += 1
            elif 0.5 < diff <= 1:
                violations["0.5-1"] += 1
            elif 1 < diff <= 1.5:
                violations["1-1.5"] += 1
            elif 1.5 < diff <= 2:
                violations["1.5-2"] += 1
            elif 2 < diff <= 2.5:
                violations["2-2.5"] += 1
            elif 2.5 < diff <= 3:
                violations["2.5-3"] += 1
            else:
                violations["3<"] += 1

    print("Total # of violations:", viol_count)
    print("Total distance of violations:", diff_sum)
    csx_func.makeNOEHist(my_path, violations)

    return viol_count

# def calcNOEviolations(args, saveShifts, my_path):
#     # parse data to restraint objects returned from pypy process
#     for data in saveShifts:
#         csx_obj.Restraint_Record(data[0], data[1], data[2], data[3],
#                                  data[4], data[5], data[6], data[7])

#     # fetch all restraint from class
#     restraints = csx_obj.Restraint_Record.all_restraints

#     PDB_coords    = csx_func.parse2dicts(args.PDB_file)
#     prev_id       = -1
#     avg_distances = {}
#     str_distaces  = {}

#     for restraint in restraints:
#         curr_id = int(restraint.curr_distID)

#         if prev_id == curr_id:
#             model_avg_dist = csx_func.getModelAvgDistance(PDB_coords,
#                                                           restraint.seq_ID1,
#                                                           restraint.atom_ID1,
#                                                           restraint.seq_ID2,
#                                                           restraint.atom_ID2)

#             avg_distances[curr_id].append(model_avg_dist)

#         else:
#             prev_id = curr_id
#             avg_distances[curr_id] = []
#             str_distaces[curr_id] = restraint.dist_max

#             model_avg_dist = csx_func.getModelAvgDistance(PDB_coords,
#                                                           restraint.seq_ID1,
#                                                           restraint.atom_ID1,
#                                                           restraint.seq_ID2,
#                                                           restraint.atom_ID2)

#             avg_distances[curr_id].append(model_avg_dist)

#     # averaging over the same restraint ID data
#     for key in avg_distances.keys():
#         avg = 0.0

#         for distance in avg_distances[key]:
#             avg += math.pow(float(distance), -6)

#         avg_distances[key] = math.pow(avg / len(avg_distances[key]), -1.0/6)

#     avg_dist_keys = avg_distances.keys()
#     avg_dist_keys.sort()
#     violations = {"0-0.5" : 0, "0.5-1" : 0, "1-1.5" : 0,
#                   "1.5-2" : 0, "2-2.5" : 0, "2.5-3" : 0, "3<" : 0}
#     viol_count = 0
#     diff_sum = 0

#     for key in avg_dist_keys:
#         if avg_distances[key] > str_distaces[key]:
#             viol_count += 1
#             diff = avg_distances[key] - str_distaces[key]
#             diff_sum += diff

#             if diff <= 0.5:
#                 violations["0-0.5"] +=1
#             elif 0.5 < diff <= 1:
#                 violations["0.5-1"] +=1
#             elif 1 < diff <= 1.5:
#                 violations["1-1.5"] +=1
#             elif 1.5 < diff <= 2:
#                 violations["1.5-2"] +=1
#             elif 2 < diff <= 2.5:
#                 violations["2-2.5"] +=1
#             elif 2.5 < diff <= 3:
#                 violations["2.5-3"] +=1
#             else:
#                 violations["3<"] +=1

#     print("Total # of violations:", viol_count)
#     print("Total distance of violations:", diff_sum)
#     csx_func.makeNOEHist(my_path, violations)


def calcNMR_Pride(pdb_models, my_path):
    # write model list text file
    pdb_models = csx_func.natural_sort(pdb_models)
    model_list = open(my_path + "model_list.txt", 'w')

    for model in pdb_models:
        #model_list.write(my_path + model + "\n")
        model_list.write(my_path + "/" + model + "\n")

    model_list.write("END\n")
    model_list.close()

    # write distance dict to text file
    restraints  = csx_obj.Restraint_Record.getPRIDE_restraints()
    pride_input = open(my_path + "pride_input.txt", 'w')

    pride_input.write("HEADER\n")

    prime_distances = restraints.keys()
    prime_distances.sort()

    for distance in prime_distances:
        pride_input.write(str(distance) + ' ' +
                          str(restraints[distance]) + '\n')

    pride_input.write("END\n")
    pride_input.close()

    # create binary database for PRIDE-NMR
    DEVNULL    = open(os.devnull, 'w')
    hhdb_log   = open(my_path + "hhdb.log", 'w')
    model_list = open(my_path + "model_list.txt", 'r')
    subprocess.call([__main__.pridedb,
                    "-D", my_path + "hhdb"  # model list
                    ],
                    stdin  = model_list,
                    stdout = DEVNULL,
                    stderr = hhdb_log)

    hhdb_log.close()
    DEVNULL.close()
    model_list.close()

    # run PRIDE-NMR
    DEVNULL      = open(os.devnull, 'w')
    pride_input  = open(my_path + "pride_input.txt", 'r')
    pride_output = open(my_path + "pride_output.txt", 'w')
    subprocess.call([__main__.pridenmr,
                     "-D", my_path + "hhdb",
                     "-d", str(56),
                     "-b", str(10),
                     "-m", str(3)],
                     stdin  = pride_input,
                     stdout = pride_output,
                     stderr = DEVNULL)

    pride_input.close()
    pride_output.close()
    DEVNULL.close()

    pride_scores = {}

    pride_output = open(my_path + "pride_output.txt", 'r')
    for line in pride_output:
        if line.startswith("PRIDENMR:"):
            model_num   = int(line.split()[-1])
            model_score = float(line.split()[1])
            pride_scores[model_num] = model_score

    scores = pride_scores.values()

    avg = sum(scores) * 1.0 / len(scores)
    variance = map(lambda x: (x - avg) ** 2, scores)
    standard_deviation = math.sqrt(sum(variance) * 1.0 / len(variance))

    PRIDE_data = []

    print("MAX: ", max(pride_scores, key=pride_scores.get))
    PRIDE_data.append(max(pride_scores, key=pride_scores.get))
    print("MIN: ", min(pride_scores, key=pride_scores.get))
    PRIDE_data.append(min(pride_scores, key=pride_scores.get))
    print("AVG: ", avg)
    PRIDE_data.append(avg)
    print("DEV: ", standard_deviation)
    PRIDE_data.append(standard_deviation)

    csx_func.makeNMRPrideGraph(my_path, scores, avg)

    return PRIDE_data
