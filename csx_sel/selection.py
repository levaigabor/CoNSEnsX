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
    def get_best_RDC_model(measure, user_sel, RDC_lists):
        model_scores = {}
        divide_by = 0.0

        for sel_data in user_sel:
            if sel_data[0] != "RDC":
                continue

            my_data = RDC_modell_data.RDC_data[sel_data[1]][sel_data[2]]
            experimental = RDC_lists[sel_data[1] - 1][sel_data[2]]

            for model_num, model in enumerate(my_data):
                if measure == "correlation":
                    calced = csx_func.calcCorrel(model, experimental)
                elif measure == "q-value":
                    calced = csx_func.calcQValue(model, experimental)
                elif measure == "rmsd":
                    calced = csx_func.calcRMSD(model, experimental)


                if model_num in model_scores.keys():
                    model_scores[model_num] += calced * sel_data[3]
                else:
                    model_scores[model_num] = calced * sel_data[3]

            divide_by += sel_data[3]

        max_value, max_loc = -1000, -1
        min_value, min_loc = 1000, -1

        print(model_scores)

        for loc, key in enumerate(model_scores.keys()):
            model_scores[key] /= divide_by

            if model_scores[key] > max_value:
                max_value = model_scores[key]
                max_loc = loc
            if model_scores[key] < min_value:
                min_value = model_scores[key]
                min_loc = loc

        if measure == "correlation":
            return max_loc
        elif measure == "q-value" or measure == "rmsd":
            return min_loc


def get_best_S2_pair(measure, S2_dict, S2_type, args):
    model_list = csx_obj.PDB_model.model_list
    scores     = {}

    for i in range(len(model_list)):
        for j in range(i + 1, len(model_list)):

            csx_obj.PDB_model.is_fitted = False

            my_models = [model_list[i], model_list[j]]

            my_S2 = csx_func.calcS2(my_models, S2_dict[S2_type], S2_type,
                                    args.fit, args.fit_range)

            if measure == "correlation":
                calced = csx_func.calcCorrel(my_S2, S2_dict[S2_type])
            elif measure == "q-value":
                calced = csx_func.calcQValue(my_S2, S2_dict[S2_type])
            elif measure == "rmsd":
                calced = csx_func.calcRMSD(my_S2, S2_dict[S2_type])

            scores[calced] = [i, j]

    if measure == "correlation":
        scores[max(scores.keys())]
    elif measure == "q-value":
        scores[min(scores.keys())]


def averageRDCs_on(models, RDC_num, RDC_type):
    """Returns a dictonary with the average RDCs for the given RDC type:
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


def averageS2_on(models, S2_dict, S2_type, args):
    """Returns a dictonary with the average S2 values for the given S2 type:
       averageS2[residue] = value"""

    model_list = csx_obj.PDB_model.model_list
    my_models  = []
    averageS2  = {}

    for model_num in models:
        my_models.append(model_list[model_num])

    csx_obj.PDB_model.is_fitted = False

    my_S2 = csx_func.calcS2(my_models, S2_dict[S2_type], S2_type,
                            args.fit, args.fit_range)

    return my_S2


def selection_on(priority, measure, pdb_models, RDC_lists, user_sel,
                 S2_dict=None, S2_type=None, args=None,
                 min_size=None, max_size=None, overdrive=None):

    if priority == "RDC":
        in_selection = [
                RDC_modell_data.get_best_RDC_model(measure, user_sel, RDC_lists)
            ]
    elif priority == "S2":
        in_selection = get_best_S2_pair(measure, S2_dict, S2_type, args)

    first_run  = True
    above_best = 0

    if measure == "correlation":
        prev_best  = -2
    else:
        prev_best = 1000

    while True:
        model_scores = {}

        for num, pdb in enumerate(pdb_models):
            if num in in_selection:
                continue

            divide_by = 0.0
            pdb_sel = [num]

            for selected in in_selection:
                pdb_sel.append(selected)

            # print("current selection: ", pdb_sel)

            for sel_data in user_sel:
                if sel_data[0] == "RDC":
                    RDC_num    = sel_data[1]
                    RDC_type   = sel_data[2]
                    RDC_weight = sel_data[3]

                    averageRDC = averageRDCs_on(pdb_sel, RDC_num, RDC_type)
                    my_RDC     = RDC_lists[RDC_num - 1][RDC_type]

                    if measure == "correlation":
                        calced = csx_func.calcCorrel(averageRDC, my_RDC)
                    elif measure == "q-value":
                        calced = csx_func.calcQValue(averageRDC, my_RDC)
                    elif measure == "rmsd":
                        calced = csx_func.calcRMSD(averageRDC, my_RDC)

                    if num in model_scores.keys():
                        model_scores[num] += calced * RDC_weight
                    else:
                        model_scores[num] = calced * RDC_weight

                    divide_by += RDC_weight

                elif sel_data[0] == "S2":
                    S2_type   = sel_data[1]
                    S2_weight = sel_data[2]

                    averageS2 = averageS2_on(pdb_sel, S2_dict, S2_type, args)
                    experimental = S2_dict[S2_type]

                    if measure == "correlation":
                        calced = csx_func.calcCorrel(averageS2, experimental)
                    elif measure == "q-value":
                        calced = csx_func.calcQValue(averageS2, experimental)
                    elif measure == "rmsd":
                        calced = csx_func.calcRMSD(averageS2, experimental)

                    if num in model_scores.keys():
                        model_scores[num] += calced * S2_weight
                    else:
                        model_scores[num] = calced * S2_weight

                    divide_by += S2_weight

        best_num = -1

        if measure == "correlation":
            best_val = -2
        else:
            best_val = 1000

        for num in model_scores.keys():
            model_score = model_scores[num] / divide_by
            if measure == "correlation" and model_score > best_val:
                best_val = model_score
                best_num = num
            elif (measure in ["q-value", "rmsd"] and model_score < best_val):
                best_val = model_score
                best_num = num

        print("prev best: " + str(prev_best) + ", current best: " + str(best_val))
        print(pdb_sel)

        # if new selection results a higher score
        if ((measure == "correlation" and best_val > prev_best) or
            (measure in ["q-value", "rmsd"] and best_val < prev_best)):

            # reset above the best threshold
            above_best     = 0
            prev_best      = best_val
            overdrive_best = -1
            in_selection.append(best_num)

            # check if selection reached the desired maximal size (if any)
            if max_size and len(in_selection) == max_size:
                print("size limit reached!")
                in_selection = [x+1 for x in in_selection]
                in_selection.sort()
                print("numbered as in PDB file:\n", in_selection)
                raise SystemExit

        # if new selection results a lower score
        else:
            # check if overdrive is enabled
            if overdrive and overdrive > above_best:
                above_best += 1
                print("\x1b[31mwe are in overdrive with \x1b[0m" + str(above_best))
                overdrive_best = best_val
                print("overdrive_best: " + str(overdrive_best))
                print("prev_best: " + str(prev_best))
                in_selection.append(best_num)

                if measure == "correlation" and overdrive_best > prev_best:
                    prev_best  = overdrive_best
                    above_best = 0
                elif measure in ["q-value", "rmsd"] and overdrive_best < prev_best:
                    prev_best  = overdrive_best
                    above_best = 0

                if overdrive == above_best:
                    if measure == "correlation" and overdrive_best < prev_best:
                        for _ in range(overdrive):
                            # print("POP")
                            del in_selection[-1]

                        print(in_selection)

                        # in_selection = [x+1 for x in in_selection]
                        # in_selection.sort()
                        # print("numbered as in PDB file:\n", in_selection)
                        raise SystemExit

                    if measure in ["q-value", "rmsd"] and overdrive_best > prev_best:
                        for _ in range(overdrive):
                            # print("POP")
                            del in_selection[-1]

                        print(in_selection)

                        # in_selection = [x+1 for x in in_selection]
                        # in_selection.sort()
                        # print("numbered as in PDB file:\n", in_selection)
                        raise SystemExit

                continue

            # check if selection reached the desired minimal size (if any)
            if min_size and len(in_selection) < min_size:
                print("we are over the peak!")
                prev_best = best_val
                in_selection.append(best_num)
                continue


            # in_selection = [x+1 for x in in_selection]
            # in_selection.sort()
            # print("numbered as in PDB file:\n", in_selection)
            raise SystemExit
