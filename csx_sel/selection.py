#!/usr/bin/env python3

from csx_libs import objects as csx_obj
from csx_libs import methods as csx_func

class RDC_modell_data(object):

    """Class for per model RDC data"""
    RDC_data = {}

    def __init__(self, RDC_list_num, RDC_type, RDC_list_data):
        if RDC_list_num in RDC_modell_data.RDC_data.keys():
            RDC_modell_data.RDC_data[RDC_list_num][RDC_type] = RDC_list_data
        else:
            RDC_modell_data.RDC_data[RDC_list_num] = {}
            RDC_modell_data.RDC_data[RDC_list_num][RDC_type] = RDC_list_data

    @staticmethod
    def get_best_model(user_sel, RDC_lists):
        # user_sel = ["RDC", 1, "0_N_H", 1], ["RDC", 2, "0_N_H", 1]

        model_scores = {}
        divide_by = 0.0

        for sel_data in user_sel:
            if sel_data[0] != "RDC":
                continue

            my_data = RDC_modell_data.RDC_data[sel_data[1]][sel_data[2]]

            for model_num, model in enumerate(my_data):
                correl  = csx_func.calcCorrel(model, RDC_lists[sel_data[1] - 1][sel_data[2]])

                if model_num in model_scores.keys():
                    model_scores[model_num] += correl * sel_data[3]
                else:
                    model_scores[model_num] = correl * sel_data[3]

            divide_by += sel_data[3]

        max_value, max_loc = 0, 0

        for loc, key in enumerate(model_scores.keys()):
            model_scores[key] /= divide_by

            if model_scores[key] > max_value:
                max_value = model_scores[key]
                max_loc = loc

        return(max_loc)


def averageRDCs_on(models, RDC_num, RDC_type):
    """Returns a dictonary with the average RDCs for a given RDC type:
       averageRDC[residue] = value"""

    averageRDC = {}

    # my_data is a list for models which contain dictonaries
    my_data = RDC_modell_data.RDC_data[RDC_num][RDC_type]

    for model_num, model in enumerate(my_data):
        if model_num not in models:
            continue

        for resnum in model:
            if resnum in averageRDC.keys():
                averageRDC[resnum] += model[resnum]
            else:
                averageRDC[resnum] = model[resnum]

    for resnum in list(averageRDC.keys()):
        averageRDC[resnum] /= len(models)

    return averageRDC


def selection_on(pdb_models, RDC_lists, user_sel, min_size=None, max_size=None):
    best_num     = RDC_modell_data.get_best_model(user_sel, RDC_lists)
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

            #print("current selection: ", pdb_sel)

            for sel_data in user_sel:
                if sel_data[0] != "RDC":
                    continue

                RDC_num    = sel_data[1]
                RDC_type   = sel_data[2]
                RDC_weight = sel_data[3]

                averageRDC = averageRDCs_on(pdb_sel, RDC_num, RDC_type)
                my_RDC = RDC_lists[RDC_num - 1][RDC_type]
                correl  = csx_func.calcCorrel(averageRDC, my_RDC)

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

            if max_size and len(in_selection) == max_size:
                print("size limit reached!")
                #in_selection = [x+1 for x in in_selection]
                in_selection.sort()
                print("numbered as in PDB file:\n", in_selection)
                raise SystemExit
        else:
            if min_size and len(in_selection) < min_size:
                print("we are over the peak!")
                print("hi, i am dori")
                prev_best = best_val
                in_selection.append(best_num)
                continue

            #in_selection = [x+1 for x in in_selection]
            in_selection.sort()
            print("numbered as in PDB file:\n", in_selection)
            raise SystemExit
