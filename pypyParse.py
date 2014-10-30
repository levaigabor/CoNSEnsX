import nmrpystar
import sys


def parseSTR(STR_file):
    star_file = open(STR_file)
    myString = ""

    for line in star_file:
        myString += line

    parsed = nmrpystar.parse(myString)

    if parsed.status != 'success':
        print('Error during STR parsing: ', parsed)
        raise SystemExit
    else:
        return parsed


my_parsed  = parseSTR(sys.argv[1])
saveShifts = my_parsed.value

saveShiftName = "CNS/XPLOR_distance_constraints_"
NOEkey        = "Gen_dist_constraint_list.Constraint_type"

for frame in saveShifts.saves.keys():
    if frame.startswith(saveShiftName):
        if NOEkey in saveShifts.saves[frame].datums.keys():
            if saveShifts.saves[frame].datums[NOEkey] == "NOE":
                loopShifts = saveShifts.saves[frame].loops[-1]
                break

restraints = []

for ix in range(len(loopShifts.rows)):   # fetch values from STR file
    row = loopShifts.getRowAsDict(ix)

    curr_distID = row["Gen_dist_constraint.ID"]
    seq_ID1     = row["Gen_dist_constraint.Seq_ID_1"]
    seq_ID2     = row["Gen_dist_constraint.Seq_ID_2"]
    atom_ID1    = row["Gen_dist_constraint.Atom_ID_1"]
    atom_ID2    = row["Gen_dist_constraint.Atom_ID_2"]
    max_val     = row["Gen_dist_constraint.Distance_upper_bound_val"]

    restraints.append([curr_distID, seq_ID1,  seq_ID2,
                                    atom_ID1, atom_ID2, max_val])

print restraints
