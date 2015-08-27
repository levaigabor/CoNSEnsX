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
