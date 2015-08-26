#!/usr/bin/env python3

from csx_libs import objects as csx_obj


def averageRDCs_on(models, RDC_num, RDC_type):
    """Returns a dictonary with the average RDCs for a given RDC type:
       averageRDC[residue] = value"""

    averageRDC = {}

    # my_data is a list for models which contain dictonaries
    my_data = csx_obj.RDC_modell_data.RDC_lists[RDC_num][RDC_type]

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
